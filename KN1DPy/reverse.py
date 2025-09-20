#reverses the order of a list at the given dimension (subscript)
def reverse (a, subscript=1):
    ndims = 1 #initially you can assume at least 1 dimension
    b = a
    while type(b[0]) == list: #if the 1st variable is also a list then a dimension is added, recurring until no longer true
        ndims+=1
        if len(b) == 0: #accounting for the case in which an empty list is encountered
            break
        b = b[0]
    if subscript > ndims:
        raise Exception('Subscript_index must be less than or equal to number of dimensions.')
    if subscript == 1: #unique case where it is reversing the 1st dim
        a=a[::-1]
        return a
    return rev_rec(a, subscript, 1)
    
#a recursive function that iterates over everything in a, and reverses everything in the specified dim
def rev_rec(a, subscript, dim_tracker):
    i = 0
    while i < len(a):
        if dim_tracker == subscript-1:
            a[i]=a[i][::-1]
        else:
            a[i] = rev_rec(a[i], subscript, dim_tracker+1)
        i+=1
    return a
