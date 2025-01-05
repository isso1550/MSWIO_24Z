print("c0|c1|c2|c3|iter_danger|sasiedztwo|pop_size|num_iter|p_rand|")
print("c0|c1|c2|c3|iter_danger|sasiedztwo|pop_size|num_iter|p_rand|\n")
i = 0
for c0 in range(10, 17, 1):
    c0 /= 10.
    for c1 in range(0, 50, 5):
        c1 /= 100.
        for c2 in range(3,8,1):
            c2 /= 10.
            for c3 in range(0, 15, 5):
                c3 /= 100.
                print(f"{c0}|{c1}|{c2}|{c3}|160|3|50|1000|0|")
                i += 1
                #print(i)