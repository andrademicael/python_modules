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
	from numpy import concatenate, arange, array
	import ipdb
	ipdb.set_trace()
	l = a.size # string size
	if l % k == 0: # n -> numer of blocks
		n = l//k
	else:
		n = (l//k) +1

	pos = array([], dtype = int)

	j = 0
	cont = 0
	# procura quais blocos tiveram numero ímpar de correções
	for i in e_pos: # realiza uma varredura nos bits corrigidos
		while True:
			if i<(j+1)*k: # faz uma contagem de quantos bits foram corrigidos de bloco (j+1)*k
				cont += 1
				break # quando um bit é pertencente ao bloco, o laço deve ser quebrado para verificação do próximo pit
			else: 
				'''
					caso o bit procurado não esteja contido no bloco (j+1)*k, duas coisas podem acontecer:
					1 - não existe nenhum bit corrigido no bloco. Logo, é necessário descobrir em que bloco o 
					próximo bit corrigido se encontra.
					2 - o número de bits corrigidos (cont) é impar e será realizada uma correção no bloco em questão
					3 - O número de bist corrigidos (cont) é par e não será realizada uma correção
				'''
				if cont == 0:
					'''
						Quando nenhum bit corrigido faz parte do bloco, é feita uma busca para encontrar 
						o bloco onde se encontra. O loop incrementa a variável j até que ela indique o bloco onde
						o bit da posição i se encontra
					'''
					while i > (j+1)*k:
						j += 1

				elif (cont % 2) == 1:
					if j < n-1:
						b[j*k:(j+1)*k], p = binary(a[j*k:(j+1)*k], b[j*k:(j+1)*k])
						pos = concatenate((pos,[j*k + p]))
					else:
						b[j*k:], p = binary(a[j*k:], b[j*k:])
						pos = concatenate((pos,[j*k + p]))
					j += 1
					cont = 0
					#break
	return [b, pos]

def recurs(a,b,pos,sig,k,step):
	'''
		realiza a correção recorrente na reconciliação do segundo passo em diante
		'siga' contem as posições originais de cada bit após as permutações
 	'''
	import ipdb
	ipdb.set_trace()

	#e_pos = siga[pos] # recebe as posições dos erros encontrados
	#e_pos.sort() # posições, na string original, dos erros corrigidos
	pos.sort()
	isig = inv_perm(sig) # isig faz a permutação inversa ao padrão em 'sig'
	#siga = sig.copy()
	for i in range(step-1,0,-1):
		a = a[isig]
		b = b[isig]
		#siga = siga[isig]
		k //= 2
	#nesse ponto, as strings foram permutadas inversamente de modo que estão na posição original
	#step = 0?
	l = a.size # string size
	if l % k == 0: # n -> numer of blocks
		n = l//k
	else:
		n = (l//k) +1
	
	b, pos = partial_dic(a,b,pos,k)
	print('A função retornou para as posições iniciais, corrigindo os erros nas posições: \n%s' % pos)

	for j in range(i,step):
		a = a[sig]
		b = b[sig]
		siga = siga[sig]
		k *= 2
		b, pos = dichotomic(a,b,k)
		# if pos.any():
		# 	b = recurs(a,b,pos,sig,siga,k,i)
	return b

def cor_var(n, m, ro, sig):
	''' 
		função que retorna duas variáveis aleatórias gaussianas correlacionadas.
		Input:
		n - comprimento dos vetores
		m - média
		ro - correlação

		Output:
		x, y - variáveis aleatórias gaussianas correlacionadas
	'''
	from numpy.random import randn
	from numpy import array, dot
	from scipy.linalg import cholesky

	xy = sig*randn(n,2) + m

	#!!! corigir !!!
	R = array([[1, ro],[ro, 1]]) #matriz de covariancia 

	L = cholesky(R)
	
	xy_cor = dot(xy,L)
	x = (xy_cor[:,0]).transpose()
	y = (xy_cor[:,1]).transpose()

	return [x, y]

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
elif __name__ == '__plt_pdf__':
	plt_pdf()
elif __name__ == '__main__':
	b2_exp()