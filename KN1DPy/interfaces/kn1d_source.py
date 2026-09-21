"""KN1D edge neutral sources: n_e particle source and associated heat sinks."""
import collections
import dataclasses
from typing import Annotated, Literal

import chex
import jax
from jax import numpy as jnp
import numpy as np
import pydantic
from torax._src import array_typing, jax_utils, state
from torax._src.config import runtime_params as runtime_params_lib
from torax._src.geometry import geometry
from torax._src.neoclassical.conductivity import base as conductivity_base
from torax._src.sources import base, source, source_profiles
from torax._src.sources import gas_puff_source as gas_puff_source_lib
from torax._src.sources import generic_ion_el_heat_source as generic_ion_el_heat_source_lib
from torax._src.sources import generic_particle_source as generic_particle_source_lib
from torax._src.sources import runtime_params as sources_runtime_params_lib
from torax._src.torax_pydantic import torax_pydantic

from KN1DPy.common import constants as CONST
from KN1DPy.kn1d_lite import KN1DLiteResults, kn1d_lite

# Hydrogenic ionization potential (13.6 eV), in Joules. Energy removed from
# the electron population per ionization event of an injected neutral.
_DEFAULT_IONIZATION_ENERGY_J = 2.182e-18

# Name under which the KN1D particle source is registered in the TORAX
# sources dict; the heat sink sources read their KN1D parameters from it.
_GAS_PUFF_SOURCE_NAME = 'gas_puff'


# pylint: disable=invalid-name
@jax.tree_util.register_dataclass
@dataclasses.dataclass(frozen=True)
class RuntimeParams(sources_runtime_params_lib.RuntimeParams):
  separatrix_neutral_density: array_typing.FloatScalar
  separatrix_neutral_energy: array_typing.FloatScalar
  # Multi-beam incident BC: energies [eV] and number-density fractions of
  # separatrix_neutral_density. Empty tuples select the single
  # separatrix_neutral_energy beam. Tuple LENGTH is static (config-fixed),
  # so `if beam_energies:` is a valid trace-time branch.
  beam_energies: tuple
  beam_fractions: tuple
  mesh_size: int
  grid_factor: float
  #ion_rate_method: str
  collisions_H2_H_EL: bool
  collisions_H_H_EL: bool
  collisions_H_P_EL: bool
  collisions_H_P_CX: bool
  simple_charge_exchange: bool


# Memoization of the KN1D solution on its plasma-profile and configuration
# inputs, so that the particle source and the heat sink sources derived from
# the same solution share a single KN1D run per timestep.
_KN1D_CACHE_MAX_SIZE = 4
_kn1d_cache: 'collections.OrderedDict[tuple, KN1DLiteResults]' = (
    collections.OrderedDict()
)


def _run_kn1d_lite_cached(
    x, mu, Ti, Te, n, vxi, incident_n0, energy_eV, mesh_size, grid_fctr, h2_h_el, h_h_el, h_p_el, h_p_cx, simple,
    fractions=None,
) -> KN1DLiteResults:
  # energy_eV may be a scalar (single beam) or a 1D array of beam energies
  # with matching number-density `fractions` of incident_n0.
  energies = np.atleast_1d(np.asarray(energy_eV, dtype=float))
  frac = (
      None if fractions is None or np.asarray(fractions).size == 0
      else np.atleast_1d(np.asarray(fractions, dtype=float))
  )
  key = (
      np.asarray(x).tobytes(),
      float(mu),
      np.asarray(Ti).tobytes(),
      np.asarray(Te).tobytes(),
      np.asarray(n).tobytes(),
      np.asarray(vxi).tobytes(),
      float(incident_n0),
      energies.tobytes(),
      b'' if frac is None else frac.tobytes(),
      int(mesh_size),
      float(grid_fctr),
      bool(h2_h_el),
      bool(h_h_el),
      bool(h_p_el),
      bool(h_p_cx),
      bool(simple),
  )
  cached = _kn1d_cache.get(key)
  if cached is not None:
    _kn1d_cache.move_to_end(key)
    return cached
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
  # x arrives center-to-edge (descending) and is reversed so KN1D's x ascends
  # from the edge, where the neutrals enter; the profiles must be reversed
  # with it or KN1D sees the core plasma at the neutral inlet.
  result = kn1d_lite(
      x=np.array(x[::-1]),
      mu=mu,
      Ti=np.array(1.0e3 * Ti)[::-1],
      Te=np.array(1.0e3 * Te)[::-1],
      n=np.array(n)[::-1],
      vxi=np.array(vxi).flatten()[::-1],
      incident_n0=float(incident_n0),
      energies_eV=list(energies),
      fractions=None if frac is None else list(frac),
      config=kn1d_config,
  )
  _kn1d_cache[key] = result
  while len(_kn1d_cache) > _KN1D_CACHE_MAX_SIZE:
    _kn1d_cache.popitem(last=False)
  return result


def calc_kn1d_lite(
    x, mu, Ti, Te, n, vxi, incident_n0, energy_eV, mesh_size, grid_fctr, h2_h_el, h_h_el, h_p_el, h_p_cx, simple,
) -> array_typing.FloatVectorCell:
  result = _run_kn1d_lite_cached(
      x, mu, Ti, Te, n, vxi, incident_n0, energy_eV, mesh_size, grid_fctr, h2_h_el, h_h_el, h_p_el, h_p_cx, simple,
  )
  return np.interp(x, result.xH, result.Sion, right=0.0)


def calc_kn1d_cx_neutral_heating(
    x, mu, Ti, Te, n, vxi, incident_n0, energy_eV, mesh_size, grid_fctr, h2_h_el, h_h_el, h_p_el, h_p_cx, simple,
) -> array_typing.FloatVectorCell:
  """Net energy transfer rate to neutrals from ions via charge exchange [W m^-3]."""
  result = _run_kn1d_lite_cached(
      x, mu, Ti, Te, n, vxi, incident_n0, energy_eV, mesh_size, grid_fctr, h2_h_el, h_h_el, h_p_el, h_p_cx, simple,
  )
  return np.interp(x, result.xH, result.EHCX, right=0.0)


def calc_kn1d_recombination_rate(
    x, mu, Ti, Te, n, vxi, incident_n0, energy_eV, mesh_size, grid_fctr, h2_h_el, h_h_el, h_p_el, h_p_cx, simple,
) -> array_typing.FloatVectorCell:
  """Volumetric recombination rate [m^-3 s^-1]."""
  result = _run_kn1d_lite_cached(
      x, mu, Ti, Te, n, vxi, incident_n0, energy_eV, mesh_size, grid_fctr, h2_h_el, h_h_el, h_p_el, h_p_cx, simple,
  )
  return np.interp(x, result.xH, result.SRecomb, right=0.0)


# KN1DLiteResults field evaluated by _kn1d_profile_callback for each of the
# point-sample functions above (kept as the callback's profile_fn keys, and
# available for direct x-space use).
_KN1D_RESULT_FIELDS = {
    calc_kn1d_lite: 'Sion',
    calc_kn1d_cx_neutral_heating: 'EHCX',
    calc_kn1d_recombination_rate: 'SRecomb',
}


def cell_average(xH, y, x_faces):
  """Conservative average of KN1D profile y(xH) over each TORAX cell.

  xH ascends from the LCFS (0); y is taken as 0 beyond xH[-1]. x_faces =
  R_out_face[-1] - R_out_face (descending, center to edge), so cell i
  spans [x_faces[i+1], x_faces[i]]. value * width reproduces KN1D's own
  integral over the cell, unlike point-sampling at the cell center.
  """
  F = np.concatenate(([0.0], np.cumsum(
      0.5 * (y[1:] + y[:-1]) * np.diff(xH))))
  Ff = np.interp(x_faces, xH, F, left=0.0, right=F[-1])
  return (Ff[:-1] - Ff[1:]) / (x_faces[:-1] - x_faces[1:])


def _kn1d_profile_callback(
    profile_fn,
    kn1d_params: 'RuntimeParams',
    geo: geometry.Geometry,
    state: state.CoreProfiles,
) -> array_typing.FloatVectorCell:
  """Evaluates a KN1D-derived profile on the TORAX cell grid via callback.

  The KN1D input x grid is measured from the LCFS (x = R_out_face[-1] -
  R_out) and the LCFS point itself (x = 0, face boundary-condition values)
  is appended, so KN1D sees the plasma-edge values at the neutral inlet;
  the earlier coupling put x = 0 at the outermost cell CENTER, half a cell
  inside the plasma edge, and never showed KN1D the boundary values.
  Arrays keep the center-to-edge order (descending x), which
  `_run_kn1d_lite_cached` reverses itself (x and profiles). The KN1D
  result field is deposited onto the TORAX cells by conservative cell
  averaging (see `cell_average`) rather than point sampling, which
  over-applied KN1D's mm-scale edge profiles wherever a cell center
  happened to sit on a peak.
  """
  field = _KN1D_RESULT_FIELDS[profile_fn]

  # Multi-beam BC when the config carries one (tuple length is static, so
  # this branch resolves at trace time); else the single-energy beam.
  if kn1d_params.beam_energies:
    energy_arg = jnp.asarray(kn1d_params.beam_energies)
    fractions_arg = jnp.asarray(kn1d_params.beam_fractions)
  else:
    energy_arg = kn1d_params.separatrix_neutral_energy
    fractions_arg = jnp.zeros((0,))

  def host(R_out, R_out_face, Ti, Te, n, vxi, Ti_b, Te_b, n_b, vxi_b,
           incident_n0, energy_eV, fractions, mesh_size, grid_fctr, h2_h_el,
           h_h_el, h_p_el, h_p_cx, simple):
    R_out, R_out_face = np.asarray(R_out), np.asarray(R_out_face)
    x = np.append(R_out_face[-1] - R_out, 0.0)
    result = _run_kn1d_lite_cached(
        x, 2.0,
        np.append(np.asarray(Ti), Ti_b),
        np.append(np.asarray(Te), Te_b),
        np.append(np.asarray(n), n_b),
        np.append(np.asarray(vxi), vxi_b),
        incident_n0, energy_eV, mesh_size, grid_fctr,
        h2_h_el, h_h_el, h_p_el, h_p_cx, simple,
        fractions=fractions,
    )
    out = cell_average(np.asarray(result.xH),
                       np.asarray(getattr(result, field)),
                       R_out_face[-1] - R_out_face)
    return out.astype(jax_utils.get_dtype())

  cell_array_shape_dtype = jax.ShapeDtypeStruct(
      shape=(geo.torax_mesh.nx,), dtype=jax_utils.get_dtype()
  )
  return jax.pure_callback(
      host,
      cell_array_shape_dtype,
      R_out=geo.R_out,
      R_out_face=geo.R_out_face,
      Ti=state.T_i.value,
      Te=state.T_e.value,
      n=state.n_e.value,
      vxi=state.toroidal_angular_velocity.value,
      Ti_b=state.T_i.face_value()[-1],
      Te_b=state.T_e.face_value()[-1],
      n_b=state.n_e.face_value()[-1],
      vxi_b=state.toroidal_angular_velocity.face_value()[-1],
      incident_n0=kn1d_params.separatrix_neutral_density,
      energy_eV=energy_arg,
      fractions=fractions_arg,
      mesh_size=kn1d_params.mesh_size,
      grid_fctr=kn1d_params.grid_factor,
      h2_h_el=kn1d_params.collisions_H2_H_EL,
      h_h_el=kn1d_params.collisions_H_H_EL,
      h_p_el=kn1d_params.collisions_H_P_EL,
      h_p_cx=kn1d_params.collisions_H_P_CX,
      simple=kn1d_params.simple_charge_exchange,
  )


def _get_kn1d_gas_puff_params(
    runtime_params: runtime_params_lib.RuntimeParams,
    requesting_source_name: str,
) -> 'RuntimeParams':
  """Fetches the KN1D parameters of the 'gas_puff' particle source."""
  gas_puff_params = runtime_params.sources.get(_GAS_PUFF_SOURCE_NAME)
  if not isinstance(gas_puff_params, RuntimeParams):
    raise ValueError(
        f"'{requesting_source_name}' requires the '{_GAS_PUFF_SOURCE_NAME}'"
        " source to be configured with the 'kn1d' model, so that both draw"
        " from the same KN1D solution."
    )
  return gas_puff_params


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
  result = _kn1d_profile_callback(calc_kn1d_lite, source_params, geo, state)
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
  # Optional multi-beam incident BC (kn1d_lite simple mode): beam energies
  # [eV] and matching number-density fractions of separatrix_neutral_density.
  # Both empty (default) selects the single separatrix_neutral_energy beam;
  # both set overrides it.
  separatrix_neutral_energies: tuple[float, ...] = ()
  separatrix_neutral_fractions: tuple[float, ...] = ()
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

  @pydantic.model_validator(mode='after')
  def _validate_beams(self) -> 'KN1DGasPuffSourceConfig':
    e = self.separatrix_neutral_energies
    f = self.separatrix_neutral_fractions
    if len(e) != len(f):
      raise ValueError(
          'separatrix_neutral_energies and separatrix_neutral_fractions must'
          f' have the same length, got {len(e)} and {len(f)}.'
      )
    if e and abs(sum(f) - 1.0) > 1.0e-6:
      raise ValueError(
          f'separatrix_neutral_fractions must sum to 1, got {sum(f)}.'
      )
    if any(x <= 0.0 for x in e) or any(x <= 0.0 for x in f):
      raise ValueError(
          'separatrix_neutral_energies and separatrix_neutral_fractions must'
          ' be positive.'
      )
    return self

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
        beam_energies=self.separatrix_neutral_energies,
        beam_fractions=self.separatrix_neutral_fractions,
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


@jax.tree_util.register_dataclass
@dataclasses.dataclass(frozen=True)
class IonizationCoolingRuntimeParams(sources_runtime_params_lib.RuntimeParams):
  ionization_energy: array_typing.FloatScalar


def calc_ionization_cooling_source(
    runtime_params: runtime_params_lib.RuntimeParams,
    geo: geometry.Geometry,
    source_name: str,
    state: state.CoreProfiles,
    calculated_source_profiles: source_profiles.SourceProfiles | None,
    unused_conductivity: conductivity_base.Conductivity | None,
) -> tuple[array_typing.FloatVectorCell, array_typing.FloatVectorCell]:
  """Electron heat sink from ionization of KN1D-sourced edge neutrals.
  """
  del geo, state  # unused
  source_params = runtime_params.sources[source_name]
  assert isinstance(source_params, IonizationCoolingRuntimeParams)
  if (
      calculated_source_profiles is None
      or 'gas_puff' not in calculated_source_profiles.n_e
  ):
    raise ValueError(
        "'kn1d_ionization_cooling' requires the 'gas_puff' ionization-rate"
        " profile to already be computed this timestep."
    )
  ionization_rate = calculated_source_profiles.n_e['gas_puff']
  electron_sink = -source_params.ionization_energy * ionization_rate
  ion_zero = jnp.zeros_like(electron_sink)
  return (ion_zero, electron_sink)


class KN1DIonizationCoolingConfig(base.SourceModelBase):
  """Electron heat sink for ionization of KN1D-sourced edge neutrals.

  Attributes:
    ionization_energy: Energy removed from the electron population per
      ionization event [J]. Defaults to 13.6 eV converted to J.
  """

  model_name: Annotated[
      Literal['kn1d_ionization_cooling'], torax_pydantic.JAX_STATIC
  ] = 'kn1d_ionization_cooling'
  ionization_energy: torax_pydantic.TimeVaryingScalar = (
      torax_pydantic.ValidatedDefault(_DEFAULT_IONIZATION_ENERGY_J)
  )
  mode: Annotated[
      sources_runtime_params_lib.Mode, torax_pydantic.JAX_STATIC
  ] = sources_runtime_params_lib.Mode.MODEL_BASED
  is_explicit: Annotated[bool, torax_pydantic.JAX_STATIC] = True

  @property
  def model_func(self) -> source.SourceProfileFunction:
    return calc_ionization_cooling_source

  def build_runtime_params(
      self,
      t: chex.Numeric,
  ) -> IonizationCoolingRuntimeParams:
    return IonizationCoolingRuntimeParams(
        prescribed_values=tuple(
            [v.get_value(t) for v in self.prescribed_values]
        ),
        mode=self.mode,
        is_explicit=self.is_explicit,
        ionization_energy=self.ionization_energy.get_value(t),
    )

  def build_source(
      self,
  ) -> generic_ion_el_heat_source_lib.GenericIonElectronHeatSource:
    return generic_ion_el_heat_source_lib.GenericIonElectronHeatSource(
        model_func=self.model_func
    )


def calc_cx_cooling_source(
    runtime_params: runtime_params_lib.RuntimeParams,
    geo: geometry.Geometry,
    source_name: str,
    state: state.CoreProfiles,
    unused_calculated_source_profiles: source_profiles.SourceProfiles | None,
    unused_conductivity: conductivity_base.Conductivity | None,
) -> tuple[array_typing.FloatVectorCell, array_typing.FloatVectorCell]:
  """Ion heat sink from charge exchange with KN1D-sourced edge neutrals.
  """
  gas_puff_params = _get_kn1d_gas_puff_params(runtime_params, source_name)
  # EHCX is the net energy transfer rate to the neutrals from charge exchange,
  # which the ion population loses (negative values mean net ion heating by
  # hot neutrals).
  cx_neutral_heating = _kn1d_profile_callback(
      calc_kn1d_cx_neutral_heating, gas_puff_params, geo, state
  )
  ion_sink = -cx_neutral_heating
  return (ion_sink, jnp.zeros_like(ion_sink))


class KN1DChargeExchangeCoolingConfig(base.SourceModelBase):
  """Ion heat sink for charge exchange with KN1D-sourced edge neutrals.

  Draws from the same KN1D solution as the 'gas_puff' particle source, which
  must be configured with the 'kn1d' model; it has no KN1D parameters of its
  own.
  """

  model_name: Annotated[
      Literal['kn1d_cx_cooling'], torax_pydantic.JAX_STATIC
  ] = 'kn1d_cx_cooling'
  mode: Annotated[
      sources_runtime_params_lib.Mode, torax_pydantic.JAX_STATIC
  ] = sources_runtime_params_lib.Mode.MODEL_BASED
  is_explicit: Annotated[bool, torax_pydantic.JAX_STATIC] = True

  @property
  def model_func(self) -> source.SourceProfileFunction:
    return calc_cx_cooling_source

  def build_runtime_params(
      self,
      t: chex.Numeric,
  ) -> sources_runtime_params_lib.RuntimeParams:
    return sources_runtime_params_lib.RuntimeParams(
        prescribed_values=tuple(
            [v.get_value(t) for v in self.prescribed_values]
        ),
        mode=self.mode,
        is_explicit=self.is_explicit,
    )

  def build_source(
      self,
  ) -> generic_ion_el_heat_source_lib.GenericIonElectronHeatSource:
    return generic_ion_el_heat_source_lib.GenericIonElectronHeatSource(
        model_func=self.model_func
    )


def calc_recombination_particle_sink(
    runtime_params: runtime_params_lib.RuntimeParams,
    geo: geometry.Geometry,
    source_name: str,
    state: state.CoreProfiles,
    unused_calculated_source_profiles: source_profiles.SourceProfiles | None,
    unused_conductivity: conductivity_base.Conductivity | None,
) -> tuple[array_typing.FloatVectorCell, ...]:
  """Particle sink for the n_e equation from recombination of KN1D edge plasma.
  """
  gas_puff_params = _get_kn1d_gas_puff_params(runtime_params, source_name)
  recombination_rate = _kn1d_profile_callback(
      calc_kn1d_recombination_rate, gas_puff_params, geo, state
  )
  return (-recombination_rate, )


class KN1DRecombinationParticleSinkConfig(base.SourceModelBase):
  """Particle sink for the n_e equation from recombination of KN1D edge plasma.

  Draws from the same KN1D solution as the 'gas_puff' particle source, which
  must be configured with the 'kn1d' model; it has no KN1D parameters of its
  own. Registers as the 'generic_particle' source.
  """

  model_name: Annotated[
      Literal['kn1d_recombination_particle_sink'], torax_pydantic.JAX_STATIC
  ] = 'kn1d_recombination_particle_sink'
  mode: Annotated[
      sources_runtime_params_lib.Mode, torax_pydantic.JAX_STATIC
  ] = sources_runtime_params_lib.Mode.MODEL_BASED
  is_explicit: Annotated[bool, torax_pydantic.JAX_STATIC] = True

  @property
  def model_func(self) -> source.SourceProfileFunction:
    return calc_recombination_particle_sink

  def build_runtime_params(
      self,
      t: chex.Numeric,
  ) -> sources_runtime_params_lib.RuntimeParams:
    return sources_runtime_params_lib.RuntimeParams(
        prescribed_values=tuple(
            [v.get_value(t) for v in self.prescribed_values]
        ),
        mode=self.mode,
        is_explicit=self.is_explicit,
    )

  def build_source(
      self,
  ) -> generic_particle_source_lib.GenericParticleSource:
    return generic_particle_source_lib.GenericParticleSource(
        model_func=self.model_func
    )


def calc_recombination_cooling_source(
    runtime_params: runtime_params_lib.RuntimeParams,
    geo: geometry.Geometry,
    source_name: str,
    state: state.CoreProfiles,
    unused_calculated_source_profiles: source_profiles.SourceProfiles | None,
    unused_conductivity: conductivity_base.Conductivity | None,
) -> tuple[array_typing.FloatVectorCell, array_typing.FloatVectorCell]:
  """Ion and electron heat sinks from recombination of KN1D edge plasma.
  """
  gas_puff_params = _get_kn1d_gas_puff_params(runtime_params, source_name)
  recombination_rate = _kn1d_profile_callback(
      calc_kn1d_recombination_rate, gas_puff_params, geo, state
  )
  # Each recombination event removes an ion/electron pair carrying the mean
  # thermal energy (3/2 kT) of its population; T_i and T_e are in keV.
  ion_sink = -1.5e3 * CONST.Q * state.T_i.value * recombination_rate
  electron_sink = -1.5e3 * CONST.Q * state.T_e.value * recombination_rate
  return (ion_sink, electron_sink)


class KN1DRecombinationCoolingConfig(base.SourceModelBase):
  """Ion and electron heat sinks for recombination of KN1D edge plasma.

  Draws from the same KN1D solution as the 'gas_puff' particle source, which
  must be configured with the 'kn1d' model; it has no KN1D parameters of its
  own.
  """

  model_name: Annotated[
      Literal['kn1d_recombination_cooling'], torax_pydantic.JAX_STATIC
  ] = 'kn1d_recombination_cooling'
  mode: Annotated[
      sources_runtime_params_lib.Mode, torax_pydantic.JAX_STATIC
  ] = sources_runtime_params_lib.Mode.MODEL_BASED
  is_explicit: Annotated[bool, torax_pydantic.JAX_STATIC] = True

  @property
  def model_func(self) -> source.SourceProfileFunction:
    return calc_recombination_cooling_source

  def build_runtime_params(
      self,
      t: chex.Numeric,
  ) -> sources_runtime_params_lib.RuntimeParams:
    return sources_runtime_params_lib.RuntimeParams(
        prescribed_values=tuple(
            [v.get_value(t) for v in self.prescribed_values]
        ),
        mode=self.mode,
        is_explicit=self.is_explicit,
    )

  def build_source(
      self,
  ) -> generic_ion_el_heat_source_lib.GenericIonElectronHeatSource:
    return generic_ion_el_heat_source_lib.GenericIonElectronHeatSource(
        model_func=self.model_func
    )
