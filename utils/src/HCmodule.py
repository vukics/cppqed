def hc(N):
    hc=[[1.], [0.,2.]]

    def ank(n, k):
        if k>n or k<0: return 0.
        elif n<len(hc): return hc[n][k]
        else: return 2.*ank(n-1,k-1)-2.*(n-1)*ank(n-2,k)

    for n in range(2,N): hc.append([ank(n,k) for k in range(n+1)])
    return hc
