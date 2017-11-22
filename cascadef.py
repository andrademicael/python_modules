def bitpar(s):
#função para calculo de paridade de uma string binária
	from numpy import sum
	
	par = sum(s) % 2

	return par

def binary(a,b):	
#função BYNARY proposta por Brassard e Salvail,
#realiza a correção de um erro entre duas strings a e b
#retorna a string b com um erro corigido e a posição do 
#erro na string	
	if bitpar(a) == bitpar(b):
		return b
		print('Os vetores apresentam paridades iguais')
	else:		
		l = a.size
		k1 = 0
		while l>1:
			if l % 2 == 0:
				if bitpar(a[k1:k1+int(l/2)]) == bitpar(b[k1:k1+int(l/2)]):
					l = int(l/2)
					k1 = k1+l
				else:
					l = int(l/2)
			else:
				if bitpar(a[k1:k1+(int(l/2)+1)]) == bitpar(b[k1:k1+(int(l/2)+1)]):
					l = int(l/2)+1
					k1 = k1+l
				else:
					l = int(l/2)+1				
		b[k1] = (b[k1] + 1) % 2

		return [b, k1]


def inv_perm(a):
#Recebe um vetor 'a' tipo inteiro que representa uma permutação qualquer
#e retorna m vetor 'ip' que realiza uma permutação inversa a 'a'

	from numpy import empty

	ip = empty(a.shape, dtype = int)
	k = 0
	for i in a:
		ip[i] = k
		k += 1
	return ip

def dichotomic(a,b,k):
# Perform a dichotomic search through the strings 'a' and 'b' for a 'k' block size.
#
# Input:
# a, b - strings beeing reconciliated. a and b must be one-dimentional arrays
# k - block size
#
# Output:
# b - string b with even parity
# pos - corrected errors positions on b
	from numpy import concatenate, array

	l = a.size # string size
	if l % k == 0: # n -> numer of blocks
		n = l//k
	else:
		n = (l//k) +1
	pos = array([], dtype = int)

	for i in range(n):
		if i < n:
			if bitpar(a[i*k:(i+1)*k]) != bitpar(b[i*k:(i+1)*k]):
				b[i*k:(i+1)*k], p = binary(a[i*k:(i+1)*k], b[i*k:(i+1)*k])
				pos = concatenate((pos,[i*k + p]))
		else:
			if bitpar(a[i*k:]) != bitpar(b[i*k:]):
				b[i*k:], p = binary(a[i*k:], b[i*k:])
				pos = concatenate((pos,[i*k + p]))
	return [b, pos]

def partial_dic(a,b,e_pos,k):
	from numpy import concatenate, arange, array

	l = a.size # string size
	if l % k == 0: # n -> numer of blocks
		n = l//k
	else:
		n = (l//k) +1
	pos = array([], dtype = int)

	for i in range(n):
		if i < n:
			if any(arange(i*k,(i+1)*k, dtype = int) == e_pos):
			
				if bitpar(a[i*k:(i+1)*k]) != bitpar(b[i*k:(i+1)*k]):
					b[i*k:(i+1)*k], p = binary(a[i*k:(i+1)*k], b[i*k:(i+1)*k])
					pos = concatenate((pos,[i*k + p]))
		else:
			if any(arange(i*k,l, dtype = int) == e_pos):
				if bitpar(a[i*k:]) != bitpar(b[i*k:]):
					b[i*k:], p = binary(a[i*k:], b[i*k:])
					pos = concatenate((pos,[i*k + p]))
	return [b, pos]

def recurs(a,b,pos,sig, siga,k,step):
	import ipdb
	ipdb.set_trace()
# realiza a correção recorrente na reconciliação do segundo passo em diante
# siga contem as posições originais de cada bit após as permutações
	e_pos = siga[pos]
	e_pos.sort()
	isig = inv_perm(sig)
	siga = sig.copy()
	for i in range(step-1,-1,-1):
		a = a[isig]
		b = b[isig]
		k //= 2
	#nesse ponto, as strings foram permutadas inversamente de modo que estão na posição original
	#step = 0
	l = a.size # string size
	if l % k == 0: # n -> numer of blocks
		n = l//k
	else:
		n = (l//k) +1
	
	b = partial_dic(a,b,e_pos,k)[0]
	i += 1

	for j in range(i,step):
		a = a[sig]
		b = b[sig]
		siga = siga[sig]
		siga = siga[sig]
		k *= 2
		b, pos = dichotomic(a,b,k)
		if pos.any():
			b = recurs(a,b,pos,sig,siga,k,i)
	return b

if __name__ == '__bitpar__':
	bitpar()
elif __name__ == '__binary__':
	binary()
elif __name__ == '__inv_perm__':
	inv_perm()
elif __name__ == '__dichotomic__':
	dichotomic()
elif __name__ == '__partial_dic__':
	partial_dic()
elif __name__ == '__recurs__':
	recurs()