
class NormalBaseFieldElement:
    __field_power = 173

    def __init__(self, vector):
        if type(vector) == str:
            vector = list(map(int, vector))
        if len(vector) == self.__field_power:
            self.vector = vector
            self.p = 2*self.__field_power + 1
        elif len(vector) < self.__field_power:
            while len(vector) != self.__field_power:
                vector.append(0)
            self.vector = vector
            self.p = 2 * self.__field_power + 1
        elif len(vector) > self.__field_power:
            self.vector = [0] * self.__field_power
            self.p = 2 * self.__field_power + 1

    def __circular_shift_right(self, i):
        shifted_vector = self.vector[self.__field_power - i:] + self.vector[:self.__field_power - i]
        return shifted_vector

    def circular_shift_left(self, i):
        shifted_vector = self.vector[i:] + self.vector[:i]
        return shifted_vector

    def __get_multiplicative_matrix(self):
        matrix = [[0] * self.__field_power for _ in range(self.__field_power)]
        for i in range(self.__field_power):
            for j in range(self.__field_power):
                if (2 ** i + 2 ** j) % self.p == 1:
                    matrix[i][j] = 1
                elif (2 ** i - 2 ** j) % self.p == 1:
                    matrix[i][j] = 1
                elif (-1 * 2 ** i + 2 ** j) % self.p == 1:
                    matrix[i][j] = 1
                elif (-1 * 2 ** i - 2 ** j) % self.p == 1:
                    matrix[i][j] = 1
                else:
                    matrix[i][j] = 0
        return matrix

    @staticmethod
    def __mult_matrix(vector, matrix):
        result = [0 for _ in range(len(matrix[0]))]
        for i in range(173):
            for j in range(173):
                result[i] += (vector[j] * matrix[j][i])
            result[i] = result[i] % 2
        return result

    @staticmethod
    def __mult_vectors(vector1, vector2):
        result = 0
        for i in range(173):
            result += vector1[i] * vector2[i]
        return result % 2

    @staticmethod
    def get_one():
        return NormalBaseFieldElement([1] * 173)

    @staticmethod
    def get_zero():
        return NormalBaseFieldElement([0] * 173)

    def add(self, another_element):
        result = [0]*self.__field_power
        for i in range(self.__field_power):
            result[i] = self.vector[i] ^ another_element.vector[i]
        return NormalBaseFieldElement(result)

    def square(self):
        result = self.__circular_shift_right(1)
        return NormalBaseFieldElement(result)

    def trace(self):
        return sum(self.vector) % 2

    def mult(self, another_element):
        multiplicative_matrix = self.__get_multiplicative_matrix()
        result = [0]*self.__field_power
        for i in range(self.__field_power):
            u = self.circular_shift_left(i)
            v = another_element.circular_shift_left(i)
            vec_dot_matrix = NormalBaseFieldElement.__mult_matrix(u, multiplicative_matrix)
            result[i] = NormalBaseFieldElement.__mult_vectors(vec_dot_matrix, v)
        return NormalBaseFieldElement(result)

    def inverse(self):
        b = self
        k = 1
        m = list(bin(self.__field_power - 1)[2:])
        for i in range(1, len(m)):
            tmp = b
            for j in range(k):
                tmp = tmp.square()
            b = b.mult(tmp)
            k *= 2
            if m[i] == '1':
                b = b.square().mult(self)
                k += 1
        return b.square()

    def power(self, n):
        result = NormalBaseFieldElement.get_one()
        for i in range(len(n.vector) - 1):
            if n.vector[i] == 1:
                result = result.mult(self)
            result = result.square()
        if n.vector[len(n.vector) - 1] == 1:
            result = result.mult(self)
        return result

    def __str__(self):
        return "".join(map(str, self.vector[:self.__field_power - list(reversed(self.vector)).index(1)]))
