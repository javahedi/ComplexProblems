




def get_spinspin_sparse(L, N, n, m, hz, sz):
    ''' 
    Create operator S_0 S_r (if sz_only=False) or S^z_n S^z_m
    if sz_only = True
    '''
    
    
    def get_site_value(state, site):
        ''' Function to get local value at a given site '''
        return (state >> site) & 1

    # Functions to create sz basis
    def first_state(N, sz):
        ''' Return first state of Hilbert space in lexicographic order '''
        n_upspins = N//2 + sz
        return (1 << n_upspins) - 1

    def next_state(state):
        '''
        Return next state of Hilbert space in lexicographic order

        This function implements is a nice trick for spin 1/2 only,
        see http://graphics.stanford.edu/~seander/bithacks.html
        # Next Bit Permutation for details
        '''
        t = (state | (state - 1)) + 1
        return t | ((((t & -t) // (state & -state)) >> 1) - 1)

    def last_state(L, sz):
        ''' Return last state of Hilbert space in lexicographic order '''
        n_upspins = N//2 + sz
        return ((1 << n_upspins) - 1) << (N - n_upspins)


    # check if sz is valid
    assert (sz <= (N // 2 + N % 2)) and (sz >= -N//2)
    
    # Create list of states with fixed sz
    basis_states = []
    state = first_state(N, sz)
    end_state = last_state(N, sz)
    while state <= end_state:
        basis_states.append(state)
        state = next_state(state)
    
    
    
    # Empty lists for sparse matrix
    rows = []
    cols = []
    data = []
    
    

    bond = (n,m)
    
    # Run through all spin configurations
    for state_index, state in enumerate(basis_states):
        
        # Apply Ising bonds
        ising_diagonal = 0
        if get_site_value(state, bond[0]) == get_site_value(state, bond[1]):
            ising_diagonal += 0/4 +hz/2
        else:
            ising_diagonal -= 0/4 +hz/2
            
        rows.append(state_index)
        cols.append(state_index)
        data.append(ising_diagonal)

        
        
       
        # Apply exchange bond
        flipmask = (1 << bond[0]) | (1 << bond[1])
        if get_site_value(state, bond[0]) != get_site_value(state, bond[1]):
             new_state = state ^ flipmask
             new_state_index = basis_states.index(new_state)
             rows.append(state_index)
             cols.append(new_state_index)
             data.append(0.25)
             
    return rows, cols, data