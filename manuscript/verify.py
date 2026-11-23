import json
from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True, eq=False)
class Gaussian:
    real: Fraction = Fraction(0)
    imag: Fraction = Fraction(0)

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Gaussian) else Gaussian(Fraction(value))

    def __add__(self, other):
        other = self.coerce(other)
        return Gaussian(self.real + other.real, self.imag + other.imag)

    __radd__ = __add__

    def __neg__(self):
        return Gaussian(-self.real, -self.imag)

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __mul__(self, other):
        other = self.coerce(other)
        return Gaussian(self.real * other.real - self.imag * other.imag,
                        self.real * other.imag + self.imag * other.real)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        norm_squared = other.real**2 + other.imag**2
        product = self * Gaussian(other.real, -other.imag)
        return Gaussian(product.real / norm_squared, product.imag / norm_squared)

    def __eq__(self, other):
        other = self.coerce(other)
        return self.real == other.real and self.imag == other.imag


def adjoint(matrix):
    return [[Gaussian(Gaussian.coerce(value).real, -Gaussian.coerce(value).imag)
             for value in column] for column in zip(*matrix)]


def check_metric(cosines, sines):
    size = len(cosines)
    eigenvalues = [Gaussian(cosine, sine) for cosine, sine in zip(cosines, sines)]
    triangular = [[eigenvalues[row] if row == column else
                   Gaussian(2 * cosines[column]) if row < column else Gaussian()
                   for column in range(size)] for row in range(size)]
    basis = [[Gaussian(Fraction(row == column)) for column in range(size)]
             for row in range(size)]
    for column in range(size):
        for row in range(column - 1, -1, -1):
            basis[row][column] = sum(triangular[row][inner] * basis[inner][column]
                                     for inner in range(row + 1, column + 1)) / (
                                         eigenvalues[column] - eigenvalues[row])
    inverse = [[Gaussian(Fraction(row == column)) for column in range(size)]
               for row in range(size)]
    for column in range(size):
        for row in range(column - 1, -1, -1):
            inverse[row][column] = -sum(basis[row][inner] * inverse[inner][column]
                                        for inner in range(row + 1, column + 1))
    assert multiply(basis, inverse) == identity(size)
    assert multiply(inverse, basis) == identity(size)
    assert multiply(triangular, basis) == multiply(basis, diagonal(eigenvalues))
    metric = multiply(adjoint(inverse), inverse)
    assert multiply(multiply(adjoint(triangular), metric), triangular) == metric
    hermitian = [[(triangular[row][column] + adjoint(triangular)[row][column]) / 2
                  for column in range(size)] for row in range(size)]
    assert hermitian == [[cosines[max(row, column)] for column in range(size)]
                         for row in range(size)]
    return {'diagonalizer': 'T S = S D', 'inverse': 'S inverse(S) = I',
            'positive_metric': 'G = inverse(S)* inverse(S)',
            'metric_identity': 'T* G T = G', 'arithmetic': 'Q(i)'}


def identity(size):
    return [[Fraction(row == column) for column in range(size)]
            for row in range(size)]


def transpose(matrix):
    return [list(column) for column in zip(*matrix)]


def multiply(left, right):
    return [[sum(left[row][inner] * right[inner][column]
                 for inner in range(len(right)))
             for column in range(len(right[0]))]
            for row in range(len(left))]


def trace(matrix):
    return sum(matrix[index][index] for index in range(len(matrix)))


def positive_semidefinite(matrix):
    if matrix != transpose(matrix):
        return False
    residual = [row[:] for row in matrix]
    for pivot in range(len(matrix)):
        value = residual[pivot][pivot]
        if value < 0:
            return False
        if value == 0:
            if any(residual[row][pivot] != 0
                   for row in range(pivot + 1, len(matrix))):
                return False
            continue
        for row in range(pivot + 1, len(matrix)):
            for column in range(pivot + 1, len(matrix)):
                residual[row][column] -= (
                    residual[row][pivot] * residual[pivot][column] / value)
    return True


def diagonal(values):
    return [[values[row] if row == column else Fraction(0)
             for column in range(len(values))]
            for row in range(len(values))]


def integer_text(value):
    digits = str(abs(value))
    if len(digits) <= 12:
        return str(value)
    terms = []
    for offset in range(0, len(digits), 6):
        end = len(digits) - offset
        start = max(0, end - 6)
        coefficient = int(digits[start:end])
        if coefficient:
            terms.append(str(coefficient) + ('*10^' + str(offset) if offset else ''))
    return ('-' if value < 0 else '') + '(' + '+'.join(reversed(terms)) + ')'


def fraction_text(value):
    return integer_text(value.numerator) + '/' + integer_text(value.denominator)


def check_matrix_examples():
    checks = []
    for size in range(1, 7):
        for singular in (False, True):
            if size == 1 and singular:
                continue
            rank = size - int(singular)
            unitary = identity(size)
            for index in range(rank - 1):
                parameter = Fraction(index + 1, 20 * size)
                cosine = (1 - parameter**2) / (1 + parameter**2)
                sine = 2 * parameter / (1 + parameter**2)
                rotation = identity(size)
                rotation[index][index] = cosine
                rotation[index + 1][index + 1] = cosine
                rotation[index][index + 1] = -sine
                rotation[index + 1][index] = sine
                unitary = multiply(rotation, unitary)
            assert multiply(transpose(unitary), unitary) == identity(size)
            modulus = diagonal([Fraction(index + 1) if index < rank else Fraction(0)
                                for index in range(size)])
            matrix = multiply(unitary, modulus)
            assert multiply(transpose(matrix), matrix) == multiply(modulus, modulus)
            hermitian = [[(matrix[row][column] + matrix[column][row]) / 2
                          for column in range(size)] for row in range(size)]
            assert positive_semidefinite(hermitian)
            assert positive_semidefinite([
                [rank * modulus[row][column] - hermitian[row][column]
                 for column in range(size)] for row in range(size)])
            factor = [[Fraction((row + 2) * (column + 1) % 5, 7)
                       for column in range(size)] for row in range(size)]
            density = multiply(transpose(factor), factor)
            for index in range(size):
                density[index][index] += Fraction(1, size + 1)
            numerator = trace(multiply(density, hermitian))
            denominator = trace(multiply(density, modulus))
            assert denominator > 0
            ratio = numerator / denominator
            assert 0 <= ratio <= rank
            lower = Fraction(1, size + 1)
            upper = lower + sum(value**2 for row in factor for value in row)
            condition_upper = upper / lower
            condition_bound = rank * condition_upper / (condition_upper + rank - 1)
            assert ratio <= condition_bound
            if size > 1:
                assert ratio < size
            assert (numerator * 37) / (denominator * 37) == ratio
            checks.append({'dimension': size, 'rank': rank,
                           'quotient': fraction_text(ratio),
                           'condition_number_enclosure_upper': fraction_text(condition_upper),
                           'quotient_bound_from_enclosure': fraction_text(condition_bound)})
    zero = [[Fraction(0)]]
    assert positive_semidefinite(zero)
    assert not positive_semidefinite([[Fraction(0), Fraction(1)],
                                     [Fraction(1), Fraction(0)]])
    assert not positive_semidefinite([[Fraction(-1)]])
    return checks


def check_repeated_singular_values():
    checks = []
    for size in range(1, 7):
        values = [Fraction(index // 2 + 1) for index in range(size)]
        count = len(set(values))
        modulus = diagonal(values)
        unitary = identity(size)
        for index in range(size - 1):
            parameter = Fraction(index + 1, 50 * size)
            cosine = (1 - parameter**2) / (1 + parameter**2)
            sine = 2 * parameter / (1 + parameter**2)
            rotation = identity(size)
            rotation[index][index] = cosine
            rotation[index + 1][index + 1] = cosine
            rotation[index][index + 1] = -sine
            rotation[index + 1][index] = sine
            unitary = multiply(rotation, unitary)
        matrix = multiply(unitary, modulus)
        assert multiply(transpose(matrix), matrix) == multiply(modulus, modulus)
        hermitian = [[(matrix[row][column] + matrix[column][row]) / 2
                      for column in range(size)] for row in range(size)]
        assert positive_semidefinite(hermitian)
        assert positive_semidefinite([[count * modulus[row][column] - hermitian[row][column]
                                       for column in range(size)] for row in range(size)])
        density = diagonal([Fraction(index + 1) for index in range(size)])
        ratio = trace(multiply(density, hermitian)) / trace(multiply(density, modulus))
        condition_bound = Fraction(count * size, size + count - 1)
        assert ratio <= condition_bound
        checks.append({'dimension': size, 'distinct_positive_singular_values': count,
                       'condition_number': size, 'quotient': fraction_text(ratio),
                       'quotient_upper_bound': fraction_text(condition_bound)})
    return checks


def check_sharpness():
    checks = []
    parameter = Fraction(1, 1000)
    delta = parameter**2
    for size in range(1, 9):
        cosines = [(1 - (index * parameter)**2) / (1 + (index * parameter)**2)
                   for index in range(1, size + 1)]
        sines = [2 * index * parameter / (1 + (index * parameter)**2)
                 for index in range(1, size + 1)]
        assert all(cosine**2 + sine**2 == 1
                   for cosine, sine in zip(cosines, sines))
        assert len(set(zip(cosines, sines))) == size
        assert all(cosines[index] > cosines[index + 1]
                   for index in range(size - 1))
        assert cosines[-1] > 0
        hermitian = [[cosines[max(row, column)] for column in range(size)]
                     for row in range(size)]
        gram = [[Fraction(0) for column in range(size)] for row in range(size)]
        extended = cosines + [Fraction(0)]
        for index in range(size):
            weight = extended[index] - extended[index + 1]
            assert weight > 0
            for row in range(index + 1):
                for column in range(index + 1):
                    gram[row][column] += weight
        assert gram == hermitian
        assert positive_semidefinite(hermitian)
        numerator = sum(map(sum, hermitian)) + delta * sum(cosines)
        denominator = size * (1 + delta)
        ratio = numerator / denominator
        assert size - Fraction(1, 100) < ratio < size
        metric = check_metric(cosines, sines)
        checks.append({'dimension': size, 't': fraction_text(parameter),
                       'delta': fraction_text(delta),
                       'complex_matrix_identities': metric,
                       'quotient_lower_bound': fraction_text(ratio),
                       'gap_below_dimension_less_than': '1/100'})
    return checks


def main():
    examples = check_matrix_examples()
    repeated = check_repeated_singular_values()
    sharpness = check_sharpness()
    print(json.dumps({'arithmetic': 'exact rational',
                      'matrix_examples': examples,
                      'repeated_singular_value_examples': repeated,
                      'sharpness_examples': sharpness,
                      'coverage': 'Finite matrix identities and inequalities; '
                                  'all-dimension statements have written proofs in note.tex.',
                      'status': 'pass'}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
