import math


def f(x):
    return math.cos(x) - 4.4 * x


def f_der(x):
    return -math.sin(x) - 4.4


def phi(x):
    return math.cos(x) / 4.4


def run_dichotomy(a, b, eps):
    n = math.ceil(math.log2((b - a) / eps))
    for _ in range(n):
        x = (b + a) / 2
        if f(a) * f(x) < 0:
            b = x
        else:
            a = x
    return n, x


def run_chord(x0, x1, m, eps):
    n = 1
    x = x1
    f0 = f(x0)
    while True:
        if abs(f(x)) / m < eps:
            break
        x = x - f(x) * (x - x0) / (f(x) - f0)
        n += 1

    return n, x


def run_moved_chord(x0, x1, m, eps):
    n = 1
    x_prev = x0
    x = x1
    while True:
        if abs(f(x)) / m < eps:
            break
        x_new = x - f(x) * (x - x_prev) / (f(x) - f(x_prev))
        x, x_prev = x_new, x
        n += 1

    return n, x


def run_newton_method(x0, m, M, eps):
    n = 0
    x = x0
    c = M / (2.0 * m)
    while True:
        prev = x
        x = x - f(x) / f_der(x)
        n += 1
        if c * (x - prev) ** 2 < eps:
            break

    return n, x


def run_modified_newton_method(x0, m, eps):
    n = 0
    x = x0
    f0 = f_der(x0)
    c = abs(f0) / m
    while True:
        prev = x
        x = x - f(x) / f0
        n += 1
        if c * abs(x - prev) < eps:
            break

    return n - 1, prev


def run_simple_iteration(a, b, q, eps):
    x0 = (b + a) / 2
    x1 = phi(x0)
    c = (1 - q) / abs(x1 - x0)
    n = math.ceil(math.log(c * eps, q))
    x = x1
    for _ in range(n - 2):
        x = phi(x)

    return n, x


def main():
    a = -math.pi / 3
    b = math.pi / 3
    eps = 0.5 * (10 ** -5)
    m = 4.4 + math.sin(a)
    M = 1
    q = math.sin(b) / 4.4
    print("Метод половиннгого деления")
    n, x = run_dichotomy(a, b, eps)
    print(f"{n:^5}|{x}")
    print("Метод неподвижных хорд деления")
    n, x = run_chord(b, a, m, eps)
    print(f"{n:^5}|{x}")
    print("Метод подвижных хорд деления")
    n, x = run_moved_chord(b, a, m, eps)
    print(f"{n:^5}|{x}")
    print("Метод Ньютона")
    n, x = run_newton_method(b, m, M, eps)
    print(f"{n:^5}|{x}")
    print("Модифицированный метод Ньютона")
    n, x = run_modified_newton_method(b, m, eps)
    print(f"{n:^5}|{x}")
    print("Метод простой итерации")
    n, x = run_simple_iteration(a, b, q, eps)
    print(f"{n:^5}|{x}")


if __name__ == "__main__":
    main()
