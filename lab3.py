import copy
import math

A = [[-0.12, 1.58, 1.63], [0.92, -0.80, 0.16], [-0.28, -0.22, -0.50]]
b = [6.23, 1.52, -2.28]
res = [2.0, 1.0, 3.0]
eps = 0.5*10e-4


#не сходится
def jacobi(a, b, eps):
    n = len(b)
    x = [0 for _ in range(0, n)]
    i = 1
    max_i = 100
    while i < max_i:
        x_prev = copy.deepcopy(x)
        for k in range(0, n):
            s = 0
            for j in range(0, n):
                if j != k:
                    s = s + a[k][j] * x[j]
            x[k] = b[k] / a[k][k] - s / a[k][k]
        if math.sqrt(sum((x[i] - x_prev[i]) ** 2 for i in range(n))) <= eps:
            break
        i += 1
    if i >= max_i:
        print('Расходится')
        return
    print('Количество итераций: ', i)
    return x


#не сходится
def seidel(A, b, eps):
    n = len(A)
    x = [0 for _ in range(0, n)]

    c = 1
    max_c = 100
    while i < max_c:
        x_new = copy.deepcopy(x)
        for i in range(n):
            s1 = sum(A[i][j] * x_new[j] for j in range(i))
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (b[i] - s1 - s2) / A[i][i]

        if math.sqrt(sum((x[i] - x_new[i]) ** 2 for i in range(n))) <= eps:
            break
        c += 1
        x = x_new
    if i >= max_c:
        print('Расходится')
        return
    print('Количество итераций: ', c)
    return x


def gauss(A, b):
    for i in range(len(A)):
        a = A[i][i]
        for k in range(i, len(A)):
            A[i][k] /= a
        b[i] /= a
        for k in range(i + 1, len(A)):
            c = -A[k][i]
            for j in range(i, len(A)):
                A[k][j] += c * A[i][j]
            b[k] += c * b[i]

    for i in range(len(A) - 1, -1, -1):
        for k in range(i - 1, -1, -1):
            c = -A[k][i]
            A[k][i] += c * A[i][i]
            b[k] += c * b[i]

    return b.copy()


def gauss_main_value(A, b):
    n = len(A)
    visited = [False] * n
    order = []
    for _ in range(n):
        max_items = map(lambda x: (x[0], *max(enumerate(x[1]), key=lambda y: abs(y[1]))), enumerate(A))
        i, j = max(filter(lambda y: not visited[y[0]], max_items), key=lambda x: abs(x[2]))[:2]
        a = A[i][j]
        visited[i] = True
        order.append((i, j))
        for k in range(n):
            A[i][k] /= a
        b[i] /= a
        for k in range(n):
            if visited[k]:
                continue
            c = -A[k][j]
            for l in range(n):
                A[k][l] += c * A[i][l]
            b[k] += c * b[i]

    x = [0 for _ in range(0, n)]
    for i in range(n - 1, -1, -1):
        for k in range(i - 1, -1, -1):
            c = -A[order[k][0]][order[i][1]]
            A[order[k][0]][order[i][1]] += c * A[order[i][0]][order[i][1]]
            b[order[k][0]] += c * b[order[i][0]]
        x[order[i][1]] = b[order[i][0]]

    return x


if __name__ == '__main__':
    a_copy = copy.deepcopy(A)
    b_copy = copy.deepcopy(b)
    print('Метод Гаусса:')
    print(gauss(a_copy, b_copy), '\n')

    a_copy = copy.deepcopy(A)
    b_copy = copy.deepcopy(b)
    print('Метод Гаусса с выбором главного элемента:')
    print(gauss_main_value(a_copy, b_copy), '\n')

    a_copy = copy.deepcopy(A)
    b_copy = copy.deepcopy(b)
    print('Метод Якоби:')
    print(jacobi(a_copy, b_copy, eps), '\n')

    a_copy = copy.deepcopy(A)
    b_copy = copy.deepcopy(b)
    print('Метод Гаусса-Зейделя:')
    print(seidel(a_copy, b_copy, eps), '\n')
