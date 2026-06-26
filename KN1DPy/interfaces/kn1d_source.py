"""KN1D edge neutral source for the n_e equation."""
import dataclasses
from typing import Annotated, ClassVar, Literal

import chex
import jax
import pydantic
import numpy as np
from torax._src import array_typing
from torax._src import jax_utils
from torax._src import state
from torax._src.config import runtime_params as runtime_params_lib
from torax._src.geometry import geometry
from torax._src.neoclassical.conductivity import base as conductivity_base
from torax._src.sources import base
from torax._src.sources import formulas
from torax._src.sources import runtime_params as sources_runtime_params_lib
from torax._src.sources import source
from torax._src.sources import source_profiles
from torax._src.sources import gas_puff_source as gas_puff_source_lib
from torax._src.torax_pydantic import torax_pydantic

from KN1DPy.kn1d_lite import kn1d_lite


# pylint: disable=invalid-name
@jax.tree_util.register_dataclass
@dataclasses.dataclass(frozen=True)
class RuntimeParams(sources_runtime_params_lib.RuntimeParams):
  separatrix_neutral_density: array_typing.FloatScalar
  separatrix_neutral_energy: array_typing.FloatScalar
  mesh_size: int
  grid_factor: float
  #ion_rate_method: str
  collisions_H2_H_EL: bool
  collisions_H_H_EL: bool
  collisions_H_P_EL: bool
  collisions_H_P_CX: bool
  simple_charge_exchange: bool


def calc_kn1d_lite(
    x, mu, Ti, Te, n, vxi, incident_n0, energy_eV, mesh_size, grid_fctr, h2_h_el, h_h_el, h_p_el, h_p_cx, simple,
) -> array_typing.FloatVectorCell:
  kn1d_config = {
      "kinetic_h": {
          "mesh_size": int(mesh_size),
          "grid_fctr": float(grid_fctr),
          "ion_rate": 'adas',
          "extra_energy_bins_eV": [],
          "ci_test": False,
          "alpha_cx_test": False,
      },
      "collisions": {
          "H2_H_EL": h2_h_el,
          "H_H_EL": h_h_el,
          "H_P_EL": h_p_el,
          "H_P_CX": h_p_cx,
          "SIMPLE_CX": simple,
      }
  }
  #out = np.zeros_like(x)
  result = kn1d_lite(
      x=np.array(x[::-1]),
      mu=mu,
      Ti=np.array(1.0e3 * Ti),
      Te=np.array(1.0e3 * Te),
      n=np.array(n),
      vxi=np.array(vxi).flatten(),
      incident_n0=float(incident_n0),
      energies_eV=[float(energy_eV)],
      config=kn1d_config,
  )
  out = np.interp(x, result.xH, result.Sion, right=0.0)
  #out[-1] = 0.0
  return out


def calc_edge_neutrals_source(
    runtime_params: runtime_params_lib.RuntimeParams,
    geo: geometry.Geometry,
    source_name: str,
    state: state.CoreProfiles,
    unused_calculated_source_profiles: source_profiles.SourceProfiles | None,
    unused_conductivity: conductivity_base.Conductivity | None,
) -> tuple[array_typing.FloatVectorCell, ...]:
  """Calculates external source term for n from puffs."""
  source_params = runtime_params.sources[source_name]
  assert isinstance(source_params, RuntimeParams)
  cell_array_shape_dtype = jax.ShapeDtypeStruct(
      shape=(geo.torax_mesh.nx,), dtype=jax_utils.get_dtype()
  )
  result = jax.pure_callback(
      calc_kn1d_lite,
      cell_array_shape_dtype,
      x=geo.R_out[-1] - geo.R_out,
      mu=2.0,
      Ti=state.T_i.value,
      Te=state.T_e.value,
      n=state.n_e.value,
      vxi=state.toroidal_angular_velocity.value,
      incident_n0=source_params.separatrix_neutral_density,
      energy_eV=source_params.separatrix_neutral_energy,
      mesh_size=source_params.mesh_size,
      grid_fctr=source_params.grid_factor,
      h2_h_el=source_params.collisions_H2_H_EL,
      h_h_el=source_params.collisions_H_H_EL,
      h_p_el=source_params.collisions_H_P_EL,
      h_p_cx=source_params.collisions_H_P_CX,
      simple=source_params.simple_charge_exchange,
  )
  return (result, )


class KN1DGasPuffSourceConfig(base.SourceModelBase):
  """Gas puff source for the n_e equation.

  Attributes:
  """

  model_name: Annotated[Literal['kn1d'], torax_pydantic.JAX_STATIC] = (
      'kn1d'
  )
  separatrix_neutral_density: torax_pydantic.TimeVaryingScalar = (
      torax_pydantic.ValidatedDefault(1.0e15)  # Density of inwardly-propagating neutrals in m**-3.
  )
  separatrix_neutral_energy: torax_pydantic.TimeVaryingScalar = (
      torax_pydantic.ValidatedDefault(3.0)     # Neutral energy at separatrix in eV.
  )
  mesh_size: pydantic.PositiveInt = 10
  grid_factor: pydantic.PositiveFloat = 0.2
  #ion_rate_method: Annotated[str, torax_pydantic.JAX_STATIC] = 'adas'
  collisions_H2_H_EL: Annotated[bool, torax_pydantic.JAX_STATIC] = True
  collisions_H_H_EL: Annotated[bool, torax_pydantic.JAX_STATIC] = True
  collisions_H_P_EL: Annotated[bool, torax_pydantic.JAX_STATIC] = True
  collisions_H_P_CX: Annotated[bool, torax_pydantic.JAX_STATIC] = True
  simple_charge_exchange: Annotated[bool, torax_pydantic.JAX_STATIC] = True
  mode: Annotated[
      sources_runtime_params_lib.Mode, torax_pydantic.JAX_STATIC
  ] = sources_runtime_params_lib.Mode.MODEL_BASED

  @property
  def model_func(self) -> source.SourceProfileFunction:
    return calc_edge_neutrals_source

  def build_runtime_params(
      self,
      t: chex.Numeric,
  ) -> RuntimeParams:
    return RuntimeParams(
        prescribed_values=tuple(
            [v.get_value(t) for v in self.prescribed_values]
        ),
        mode=self.mode,
        is_explicit=True,
        separatrix_neutral_density=self.separatrix_neutral_density.get_value(t),
        separatrix_neutral_energy=self.separatrix_neutral_energy.get_value(t),
        mesh_size=self.mesh_size,
        grid_factor=self.grid_factor,
        #ion_rate_method=self.ion_rate_method,
        collisions_H2_H_EL=self.collisions_H2_H_EL,
        collisions_H_H_EL=self.collisions_H_H_EL,
        collisions_H_P_EL=self.collisions_H_P_EL,
        collisions_H_P_CX=self.collisions_H_P_CX,
        simple_charge_exchange=self.simple_charge_exchange,
    )

  def build_source(self) -> gas_puff_source_lib.GasPuffSource:
    return gas_puff_source_lib.GasPuffSource(model_func=self.model_func)
