print("accept_rate|pa_rate|pa_band_coef|pop_size|num_iter|p_rand|")
print("accept_rate|pa_rate|pa_band_coef|pop_size|num_iter|p_rand|\n")
i = 0
#for rastrigin 10
for racc in range(60, 90+10, 5):
    racc /= 100.
    for rpa in range(50, 90+10, 5):
        rpa /= 100.
        print(f"{racc}|{rpa}|200|30|500|0")
        i += 1
        #print(i)