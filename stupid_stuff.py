import re
from itertools import chain

text = '7,21,22,28,33,37-39,52,59,62,63,65-71,74,75,77-84,88-91,93-100,102-104,106,108,110,113,119,122,123,132-137,139,140,144'
numbers = re.split(r'[,\-]', text)
seperators = re.split(r'\d+', text)[1:]
new_numbers = [str(int(num)-1) for num in numbers]
complete_list = list(chain.from_iterable(zip(new_numbers, seperators)))
print(''.join(complete_list))