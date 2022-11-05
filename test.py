from time import sleep
x = 10
for i in range(x):
    print(f'Progress {i/x:4.0%}', end='\r')
    sleep(1)