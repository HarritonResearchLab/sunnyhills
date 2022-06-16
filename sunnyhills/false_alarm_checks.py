
### ACTUAL TESTS BELOW ###

def tls_even_odd(tls_results): 
    eb_flag = True 

    odd = tls_results.depth_mean_even 

    odd = [odd[0]-odd[1], odd[0]+odd[1]]

    even = tls_results.depth_mean_odd
    even = [even[0]-even[1], even[0]+even[1]]

    temp = np.sort([even, odd])

    if temp[0][1]>temp[1][0]: 
        eb_flag = False

    print(temp)

    return eb_flag
