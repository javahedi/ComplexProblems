
def get_site_value(state, site):
    ''' Function to get local value at a given site '''
    return (state >> site) & 1

def hilbertspace_dimension(N):
    ''' return dimension of hilbertspace '''
    return 2**N



def get_spinspin_sparse(L, N, n, m, hz, sz_only=True):
    ''' 
    Create operator S_0 S_r (if sz_only=False) or S^z_n S^z_m
    if sz_only = True
    '''
    # Empty lists for sparse matrix
    rows = []
    cols = []
    data = []

    bond = (n,m)
    
    # Run through all spin configurations
    for state in range(hilbertspace_dimension(N)):
        
        # Apply Ising bonds
        ising_diagonal = 0
        if get_site_value(state, bond[0]) == get_site_value(state, bond[1]):
            ising_diagonal += 0/4 +hz/2
        else:
            ising_diagonal -= 0/4 +hz/2
        rows.append(state)
        cols.append(state)
        data.append(ising_diagonal)

        
        if not sz_only:
            # Apply exchange bond
            flipmask = (1 << bond[0]) | (1 << bond[1])
            if get_site_value(state, bond[0]) != get_site_value(state, bond[1]):
                new_state = state ^ flipmask
                rows.append(new_state)
                cols.append(state)
                data.append(0.25)

    return rows, cols, data