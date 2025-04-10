#%% Gerando vetor p
import numpy as np
import random


def generate_vector_p():
    p = [0.5, 0.2, 0.1]
    q = [0.5, 0.2, 0.1]
    for i in range(1,5):
      for j in range(3):
        q[j] = q[j]/10
        p.append(q[j])
    return p

vector_p = generate_vector_p()
print(vector_p)


#%% Criando os bits para hamming
def generate_random_bits(n):
  return np.random.randint(0, 2, n, dtype=np.uint8)

random_bits = generate_random_bits(100)
random_bits

def split_bits_reshape(total_bits, k):
  return total_bits.reshape(-1, k)

split_arrays = split_bits_reshape(random_bits, 4)
# print(split_arrays[:5])

#%% Codificador de Hamming
def hamming_codificator(split_arrays):
    encoded_arrays = []

    for arr in split_arrays:
        arr = arr.tolist()
        p1 = (arr[0] + arr[1] + arr[2]) % 2
        p2 = (arr[0] + arr[2] + arr[3]) % 2
        p3 = (arr[0] + arr[1] + arr[3]) % 2

        arr.append(p1)
        arr.append(p2)
        arr.append(p3)

        encoded_arrays.append(arr)

    return np.array(encoded_arrays)

encoded_data = hamming_codificator(split_arrays)
# print(encoded_data)

#%% Canal BSC
def BSC(encoded_data, p):
    modified_matrix = []
    for word in encoded_data:
        modified_word = [(1 + x) % 2 if random.random() < p else x for x in word]
        modified_matrix.append(modified_word)
    return np.array(modified_matrix)

p_value = 0.1
modified_matrix = BSC(encoded_data, p_value)

