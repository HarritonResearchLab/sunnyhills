def template_function(arg_one: int, arg_two: int = 5): 
    
    '''
    Args:
        arg_one: a integer
        arg_two: another integer, default is 5
    Returns: 
        arr: numpy array made from [arg_one, arg_two] 
    '''

    # imports 
    import numpy as np

    # code
    arr = np.array([arg_one, arg_two])
    

    return arr
