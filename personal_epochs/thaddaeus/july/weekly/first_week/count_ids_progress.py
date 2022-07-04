
def message(time_so_far): 
    counter = 0
    with open('personal_epochs/thaddaeus/july/weekly/first_week/log.txt', 'r') as f: 
        for line in f: 
            if '[' not in line: 
                counter+=1

    percent_done = 100*counter/10215
    etd = (10215-counter)*time_so_far/counter 
    print('finished: ', counter, 'percent complete: ', round(percent_done, 1), 'ETD: ', round(etd), round(etd/60, 2), '(minutes, hours)')

message(139)
#finished:  199 percent complete:  1.9 ETD:  617 10.28 (minutes, hours)