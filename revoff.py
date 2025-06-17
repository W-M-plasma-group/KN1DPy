def revoff():
    b = [0] * 4
    b[0] = bytes('1B')
    b[1] = bytes('[')
    b[2] = bytes('0')
    b[3] = bytes('m')
    b = str(b)
    print(b[:3])
