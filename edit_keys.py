def edit_keys(dict):
    # changes the keys for the C-Mod input 
    # Inputs:
    #   dict - the dictionary outputted from the C-Mod Save file
    # Ouputs:
    #   dict - the ammended dictionary 
    # Gwendolyn Galleher 
    dict['xlimiter'] = dict['x_lim']
    del dict['x_lim']
    dict['xsep'] = dict['x_sep']
    del dict['x_sep']
    dict['Ti'] = dict['t_i']
    del dict['t_i']
    dict['Te'] = dict['t_e']
    del dict['t_e']
    dict['n'] = dict['n_e'] #  this could be wrong 
    del dict['n_e']
    dict['vxi'] = dict['vx']
    del dict['vx']
    dict['LC'] = dict['lc']
    del dict['lc']
    dict['PipeDia'] = dict['d_pipe']
    del dict['d_pipe']
    return dict

