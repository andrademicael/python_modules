'''
	Módulo com funções para implementação do protocolo cascade e geração de duas 
	sequencias de bits correlacionadas, a partir expanção em base 2 da função 
	Fx(x) de uma pdf gaussiana
'''
def bitpar(s):
	'''
		função para calculo de paridade de uma string binária
	'''
	from numpy import sum
	par = sum(s) % 2
	return par

def binary(a,b):
	'''
		função BYNARY proposta por Brassard e Salvail,
		realiza a correção de um erro entre duas strings a e b
		retorna a string b com um erro corigido e a posição do 
		erro na string	
	'''
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
	'''
		Recebe um vetor 'a' tipo inteiro que representa uma permutação qualquer
		e retorna m vetor 'ip' que realiza uma permutação inversa a 'a'
	'''
	from numpy import empty

	ip = empty(a.shape, dtype = int)
	k = 0
	for i in a:
		ip[i] = k
		k += 1
	return ip

def dichotomic(a,b,k):
	'''
		Perform a dichotomic search through the strings 'a' and 'b' for a 'k' block size.

		Input:
		a, b - strings beeing reconciliated. a and b must be one-dimentional arrays
		k - block size

		Output:
		b - string b with even parity
		pos - corrected errors positions on b
	'''
	from numpy import concatenate, array

	l = a.size # string size
	if l % k == 0: # n -> numer of blocks
		n = l//k
	else:
		n = (l//k) +1
	pos = array([], dtype = int)

	for i in range(n):
		if i < n-1:
			if bitpar(a[i*k:(i+1)*k]) != bitpar(b[i*k:(i+1)*k]):
				b[i*k:(i+1)*k], p = binary(a[i*k:(i+1)*k], b[i*k:(i+1)*k])
				pos = concatenate((pos,[i*k + p]))
		else:
			if bitpar(a[i*k:]) != bitpar(b[i*k:]):
				b[i*k:], p = binary(a[i*k:], b[i*k:])
				pos = concatenate((pos,[i*k + p]))
	return [b, pos]

def partial_dic(a,b,e_pos,k):
	'''
		Função utilizada para recursividade no CASCADE. 
		A função está instável e, conforme conversa com o professor Bruno, não será utilizada recursividade.
	'''
	from numpy import concatenate, arange, array, zeros

	l = a.size # string size
	if l % k == 0: # n -> numer of blocks
		n = l//k
	else:
		n = (l//k) +1

	pos = array([], dtype = int)

	teste = zeros((1, l), dtype = int)
	teste[e_pos] = 1

	for j in range(n):
		if j < n-1:
			if bitpar(teste[j*k:(j+1)*k]) == 1:
				b[j*k:(j+1)*k], p = binary(a[j*k:(j+1)*k], b[j*k:(j+1)*k])
				pos = concatenate((pos,[j*k + p]))
		else:
			if bitpar(teste[j*k:]) == 1:
				b[j*k:], p = binary(a[j*k:], b[j*k:])
				pos = concatenate((pos,[j*k + p]))
	return [b, pos]

def recurs(a,b,pos,sig,k,step):
	'''
		realiza a correção recursiva na reconciliação, do segundo passo em diante.
		'siga' contém as posições originais de cada bit após as permutações
		'recurs' sempre restaurará as strings para suas condições iniciais. Logo,
		'pos' deve informar a posição original do bit corrigido.
			
		A função está instável e, conforme conversa com o professor Bruno, não será utilizada recursividade.
 	'''
	from numpy import arange
	import ipdb
	ipdb.set_trace()

	pos.sort()
	isig = inv_perm(sig) # isig faz a permutação inversa ao padrão em 'sig'
	for i in range(step-1,0,-1):
		a = a[isig]
		b = b[isig]
		k //= 2
	#nesse ponto, as strings foram permutadas inversamente de modo que estão na posição original. i = 1.

	l = a.size # string size
	siga = arange(l)
	if l % k == 0: # n -> numer of blocks
		n = l//k
	else:
		n = (l//k) +1
	
	b, pos = partial_dic(a,b,pos,k) # essa operação é semelhate à 'dichotomic' na linha 78, 

	k *= 2
	for j in range(i,step):
		a = a[sig]
		b = b[sig]
		siga = siga[sig]
		b, pos = dichotomic(a,b,k)
		if pos.any():
			b = recurs(a,b,siga[pos],sig,k,j)
		k *= 2
	return b

def cor_var(n, ro):
	''' 
		função que retorna duas variáveis aleatórias gaussianas correlacionadas.
		Input:
		n - comprimento dos vetores
		m - média
		ro - correlação

		Output:
		x, y - variáveis aleatórias gaussianas correlacionadas
		sig - variância do ruído
	'''
	from numpy.random import randn
	from numpy import array, dot, var
	from scipy.linalg import cholesky

	xy = randn(n,2)

	#!!! corigir !!!
	R = array([[1, ro],[ro, 1]]) #matriz de covariancia 

	L = cholesky(R)
	
	xy_cor = dot(xy,L)
	x = (xy_cor[:,0]).transpose()
	y = (xy_cor[:,1]).transpose()
	n = x-y

	var_x = var(x)
	var_y = var(y)
	var_n = var(n)

	return [x, y, var_x, var_y, var_n]

def strings(ro, nbits, nr):
	'''
		Duas variaveis aleatŕoia são geradas: Gaussianas correlacionadas. As funções cumulativas de 
		probabilidade de suas realizações assumem uma distribuição uniforme (maximizando entropia).
		A partir disto, são geradas expansões em base 2, dado um numero fixo de bits para representação
		dos valores das CDF's
	'''

	from numpy import log10, zeros, where
	from scipy.stats import norm

	A, B, var_a, var_b, var_n = cor_var(nr, ro) # gera as duas VA's Gaussianas correlacionadas
	
	k = int(0.73/var_n) # First iteration block size
	snr = 10*log10(var_a/var_n)

	cdfa = norm.cdf(A) # valor da CDF para cada valor de A
	cdfb = norm.cdf(B) # valor da CDF para cada valor de B

	a = zeros((nr, nbits), dtype = int)
	b = zeros((nr, nbits), dtype = int)
	for i in range(nr):
		a[i,:] = b2_exp(cdfa[i],nbits)
		b[i,:] = b2_exp(cdfb[i],nbits)

	a = a.reshape((1,a.size))[0]
	b = b.reshape((1,b.size))[0]

	print('Parâmetros da simulação: ')
	print('\nCoeficiente de correlação das VA\'s: %.2f' % ro)
	print('Realizações: %i ' % nr)
	print('Número de bits para representação: %i' % nbits)
	print('Comprimento total das strings: %i' % a.size)
	print('Quantidade de erros gerados: %i ' % where(a != b)[0].size)
	print('Variância do ruído: %.2f' % var_n)
	print('SRN: %.2f' % snr)
	print('Tamanho to bloco: %i' % k)
	
	return [a, b, k, snr]

def plt_pdf(x, *args, **kargs):
	from numpy import linspace
	from scipy.stats.kde import gaussian_kde

	kdex = gaussian_kde(x)
	dist_space = linspace( min(x), max(x), 500 )

	pdfx = kdex(dist_space)
	
	if 'plot' in kargs:
		from matplotlib.pyplot import plot, show, legend

		plot(dist_space, pdfx)
		legend(['X'])
		show()

	return [pdfx]

def b2_exp(n,nb):
	'''
		Recebe um número real entre no intervalo (0,1) e a precisão da expansão
		e retorna um string com os bits referentes à representação
	'''

	from numpy import array, concatenate, sign
	expanse = array([], dtype = int)
	flag = 0
	i = 1

	while flag == 0:
		k = 2**-i
		if sign(n-k)+1:
			b = [1]
			n -= k
		else:
			b = [0]
		expanse = concatenate((expanse,b))
		i +=1
		if expanse.size == nb:
			flag = 1
	return expanse

def cascade(a, b, sig, k):
	from numpy import where, arange

	print('\nInício da CASCADE')

	######################################################################################################
	#									CASCADE 1st Step 			           		                     #
	#                                   													             #

	# first strings permutation
	a = a[sig]
	b = b[sig] # a primeira permutação é considerada como posição inicial da cadeia de bits
	siga = arange(a.size) # original bits positions

	b, pos = dichotomic(a,b,k) # first search
	'''
		pos inform wrong bits positions that were solved on a not permuted list, e.g, arange(1,l)
		the original solved bits positions rely on siga[pos]
	'''
	step = 1 # first step compleated
	print('')
	print('Após o passo %i, restam %i erros a serem corrigidos.' % (step, where(a != b)[0].size))

	######################################################################################################
	#									CASCADE 2nd Step 			           		                     #
	#                                   													             #
	# k *= 2 # double block size
	while step<10 and where(a != b)[0].size != 0: # o tamanho do bloco para realizar dichotomic() deve ser menor que o tamanho da string
		a = a[sig] # permutation
		b = b[sig] # permutation
		siga = siga[sig] # keeps the original bits positions by permuting siga
		b, pos = dichotomic(a,b,k) # realize dichotomic() 
		step += 1 # step compleated
		print('')
		print('Após o passo %i, restam %i erros a serem corrigidos.' % (step, where(a != b)[0].size))
		'''
			if pos.size != 0: # if any error was corrected, must be done a recursive search
				print('Strings antes da recursividade: \nA: \n%s \nB:\n%s' % (a.reshape((n_r, nbits)), b.reshape((n_r, nbits))))
				b = recurs(a, b, siga[pos], sig, k, step)
				print('Strings depois da recursividade: \nA: \n%s \nB:\n%s' % (a.reshape((n_r, nbits)), b.reshape((n_r, nbits))))
			k *= 2 # double block size
		'''
	return [a, b]

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
elif __name__ == '__cor_var__':
	cor_var()
elif __name__ == '__strings__':
	strings()
elif __name__ == '__plt_pdf__':
	plt_pdf()
elif __name__ == '__main__':
	b2_exp()
elif __name__ == '__cascade__':
	cascade()