with open ('rosalind_fib.txt') as f:
    n,k = map(int, f.read().strip().split())
print(n,k)
def rabbits(n,k):
    a,b = 1, 1
    for i in range(n - 2):
        a, b = b, b + k * a
    return b
print(rabbits(n,k))
def fib(n,k):
    if n == 1 or n == 2:
        return 1
    else:
        return fib(n-1, k) + k * fib(n-2, k)
print(fib(n,k))