# https://linus-mk.hatenablog.com/entry/2020/02/23/225258の記事参考（というかほぼそのまま）


# mod対応組合せ関数 (計算量: O(r))
def mod_comb(n, r, mod=998244353):
    if n < r or r < 0:
        return 0
    if n-r < r:
        r = n-r

    res = 1
    for i in range(r):
        res *= n-i
        res *= pow(i+1, -1, mod)
        res %= mod

    return res
