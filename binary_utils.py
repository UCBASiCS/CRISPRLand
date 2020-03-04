from __future__ import division
from __future__ import print_function

from scipy.special import comb

from numpy.random import choice

def binary_inprod(a,b):
    """
    Counts the (pseudo) inner product over Galois F2.

    Args:
        a: Integer 1.
        b: Integer 2.

    Returns:
        Inner product.
    """
    v = a & b
    num_ones = count_ones_in_int(v)
    return(num_ones % 2)


def bin_to_dec(s):
    return( int(s[2:],2) )


def bit_rotate_left(val,r_bits,max_bits):
    return (val << r_bits%max_bits) & (2**max_bits-1) | \
    ((val & (2**max_bits-1)) >> (max_bits-(r_bits%max_bits)))


def count_ones_in_int(x):
    """
    Counts the number of ones in an integer.

    Args:
        x: Integer.

    Returns:
        Number of ones in an integer.
    """
    return bin(x).count("1")


def local_to_global(a,u):
    # a is the local index
    # u is the map
    
    # we need a <= u
    num_ones_in_u = bin(u).count("1")
    assert a < (2**num_ones_in_u)
    
    u_binary = bin(u)
    
    
    l = len(u_binary) - 1
    r = 0

    shift = 0
    while u_binary[l] != 'b':
        if u_binary[l] == '1':
            q = a & 1
            a = (a >> 1)
            r += (q << shift)
        shift += 1
        l -= 1

    return(r)

def global_to_local(a,u):
    # a is the global index
    # u is the map
    
    u_binary = bin(u)
    
    l = len(u_binary) - 1
    
    r = 0
    shift = 0
    
    while u_binary[l] != 'b':
        q = a & 1
        a = a >> 1
        
        if u_binary[l] == '1':    
            r += q << shift
            shift += 1
        l -= 1
    
    return(r)


# generate all bit strings having d ones
def next_string_with_same_num_ones(v):
    t = (v | (v-1))+ 1
    w = t | ((( (t & -t) // (v & -v) ) >> 1) - 1 )
    return w


def all_strings_with_k_ones(bit_length,k):
    num_total = int( comb(bit_length,k) )
    c = 2**k - 1
    my_list = []
    for i in range(num_total):
        my_list.append(c)
        if i != num_total - 1:
            c = next_string_with_same_num_ones(c)
        
    return my_list

def all_strings_up_to_k_ones(bit_length,k):
    my_list = []
    
    for i in range(k+1):
        my_list = my_list + all_strings_with_k_ones(bit_length,i)
        
    return my_list


def get_random_upto_k_degree(n,k):
    # number of combinations having i ones in length n
    p = [comb(n,i) for i in range(0,k+1)]
    # normalization constant to turn into probabilities
    z = sum(p)
    # probability of having i ones in n
    choice_probabilities = [i/z for i in p]
    
    # first choose a degree
    chosen_degree = choice(k+1,1,p=choice_probabilities)[0]
    
    r = get_random_with_k_ones_over_n(n,chosen_degree)
    
    return r

def get_random_with_k_ones_over_n(n,k):
    # degree is k
    # num locations is n
    r = 0
    
    # then choose a distribution of ones depending on the degree
    # if degree is non-zero then get the ones in place
    if k != 0:
        # chose indices to set to 1
        chosen_indices = choice(n,k,replace=False)
        # turn it into integer
        for t in chosen_indices:
            r += (2**t)
    return r


def binary_list_to_dec_list(A):
    # example input ['0b110000','0b001100']
    return map(bin_to_dec,A)    


def aliased_bin(i,aliasing_patterns):
    # both inputs need to be decimal    
    y = map(lambda u: binary_inprod(i,u),aliasing_patterns)
    
    resulting_bin = 0
    r = 1
    for t in y:
        resulting_bin += t*r
        r = r << 1
        
    return resulting_bin



