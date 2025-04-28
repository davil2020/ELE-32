import numpy as np
import random

def generate_vector_p():
    p = [0.5, 0.2, 0.1]
    q = [0.5, 0.2, 0.1]
    for _ in range(1, 5):
        for j in range(3):
            q[j] = q[j] / 10
            p.append(q[j])
    return p

def generate_random_bits(num_bits):
    return np.random.randint(0, 2, num_bits)

def split_bits_reshape(total_bits, k):
    return total_bits.reshape(-1, k)

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

def BSC(encoded_data, p):
    modified_matrix = []
    for word in encoded_data:
        modified_word = [(1 + x) % 2 if random.random() < p else x for x in word]
        modified_matrix.append(modified_word)
    return np.array(modified_matrix)

def matrix_multiplication(transmited_data, Ht):
    sindromes = []
    for data in transmited_data:
        result = np.dot(data, Ht) % 2
        sindromes.append(result)
    return np.array(sindromes)

def syndrome_to_bit_position(syndrome, Ht):
    if np.all(syndrome == 0):
        return None
    for pos in range(Ht.shape[0]):
        if np.array_equal(syndrome, Ht[pos]):
            return pos
    return None

def correct_errors(word, bit_position):
    corrected_word = word.copy()
    if bit_position is not None:
        corrected_word[bit_position] = (corrected_word[bit_position] + 1) % 2
    return corrected_word

def correct_matrix(modified_matrix, Ht):
    corrected_matrix = []
    syndromes = matrix_multiplication(modified_matrix, Ht)
    for i, word in enumerate(modified_matrix):
        bit_pos = syndrome_to_bit_position(syndromes[i], Ht)
        corrected_word = correct_errors(word, bit_pos)
        corrected_matrix.append(corrected_word)
    return np.array(corrected_matrix)

def report_differences(encoded_data, corrected_matrix):
    divergence = (encoded_data + corrected_matrix) % 2
    info_bits = divergence[:, :4]
    total_errors = np.sum(info_bits)
    total_bits = info_bits.size
    error_probability = np.float64(total_errors) / np.float64(total_bits)
    return error_probability

def hamming_simulation():
    Ht = np.array([
        [1, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
        [0, 1, 1],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ])

    vector_p = generate_vector_p()
    erro_prob = []
    
    for p in vector_p:
        random_bits = generate_random_bits(1000000)
        split_arrays = split_bits_reshape(random_bits, 4)
        encoded_data = hamming_codificator(split_arrays)
        modified_matrix = BSC(encoded_data, p)
        corrected_matrix = correct_matrix(modified_matrix, Ht)
        erro = report_differences(encoded_data, corrected_matrix)
        erro_prob.append(erro)
        print(f"p={p:.5f}, BER={erro:.5e}")

    return erro_prob
