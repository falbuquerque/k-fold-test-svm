###############################################################################
# Author: Felipe Albuquerque
###############################################################################

#
# Criação da classe que representa um k-grupo. Guarda uma matriz de dados no atributo "grupo"
#
setClass("k.grupo", representation(grupo = "matrix"))

#
# Valida o valor de k, que deve estar entre 1 e N. 
# Retorna 1 quando k é inválido e 0 quando é válido
#
isKValido = function(k, n) {((k > 1) & (k < n))}

#
# Divide os dados em k grupos
# 
# dados: os dados que devem ser divididos
# k: o número de grupos. Deve estar entre 1 e o número de dados
#
dividirEmKGrupos = function(dados, k) {
    
    # Validação de k
    linhas = length(dados[,1])
    
    if (!isKValido(k, linhas)) {
    	stop(paste("O valor de k deve estar entre 1 e", linhas))
    }
    
    # Inicialização do retorno. Criação das estruturas dos grupos
    cols = length(dados[1,])
    conjuntoRetornado = c()
    
    for (j in 1:k) {
        conjuntoRetornado = c(conjuntoRetornado, new("k.grupo", grupo = matrix(nrow = 0, ncol = cols)))
    }

    # Criação dos grupos    
    
    for (i in 1:linhas) {
        pos = i %% k # Seleção do grupo baseada no resto da divisão
        
        if (pos == 0) { # Como o índice começa em "1", utilizar "k" ao invés de "0"
            pos = k
        }
        
        # Adiciona o dado ao grupo selecionado
        conjuntoRetornado[[pos]]@grupo = rbind(conjuntoRetornado[[pos]]@grupo, dados[i,])
    }
    
    conjuntoRetornado
}

#
# Retorna a matriz sem a coluna de classe
# 
# m: a matriz
# indiceColunaClasse: o índice da coluna de classe
#
extrairMatrizSemColunaDeClasse = function(m, indiceColunaClasse) {
    linhas = length(m[,1])
    cols = length(m[1,])
    
    if (indiceColunaClasse == cols) { # coluna de classe é a última
        m[1:linhas, (1:cols - 1)]
    } else {
        x = m[1:linhas, 1:(indiceColunaClasse -1)]
        cbind(m[1:linhas, (indiceColunaClasse + 1):cols])
    }
    
}

#
# Retorna a coluna de classe da matriz
# 
# m: a matriz
# indiceColunaClasse: o índice da coluna de classe
#
extrairColunaDeClasse = function(m, indiceColunaClasse) {
    linhas = length(m[,1])
    as.matrix(m[1:linhas, indiceColunaClasse:indiceColunaClasse], nrow = linhas, ncol = 1, byrow = TRUE)
}

#
# Avalia o classificador baseado em SVM utilizando a técnica de cross-validation
# 
# dados: os dados que devem ser divididos
# k: o número de grupos. Deve estar entre 1 e o número de dados
# kernel: a função kernel. O default é kernel = "radial"
# c: a taxa de erros aceita. O default é c = 1
#
avaliarSVM = function(dados, k, indiceColunaClasse, kernel = "radial", taxaErro = 1) {
    m = as.matrix(dados)
    linhas = length(m[,1])
    cols = length(m[1,])
    
    # Validação de k
    if (!isKValido(k, linhas)) {
        stop(paste("O valor de k deve estar entre 1 e", linhas))
    }

    # Execução da avaliação
    grupos = dividirEmKGrupos(m, k)
    
    print(paste("Validação cruzada para k =", k))
    print(paste("Função kernel escolhida:", kernel))
    print(paste("Valor de C:", taxaErro))
    
    tp = 0
    fp = 0
    fn = 0
    tn = 0
    total = 0
    taxaFalsaAceitacao = 0
    taxaFalsaRejeicao = 0
    sensibilidade = 0
    especificidade = 0
    erroMedio = 0
    
    for (i in 1:k) {
        
        # Escolhe o grupo de testes
        grupoTeste = i
        conjuntoTestes = grupos[[grupoTeste]]@grupo
        
        # Monta o conjunto de treinamento com k - 1 elementos
        conjuntoTreinamento = matrix(nrow = 0, ncol = cols)
        
        for (j in 1:k) {
            
            if (j != grupoTeste) {
                conjuntoTreinamento = rbind(conjuntoTreinamento, grupos[[j]]@grupo)
            }
            
        }
        
        # Criação do modelo com o conjunto de treinamento
        x = extrairMatrizSemColunaDeClasse(conjuntoTreinamento, indiceColunaClasse)
        y = extrairColunaDeClasse(conjuntoTreinamento, indiceColunaClasse)
        modelo = svm(x = x, y = y, kernel = kernel, type = "C-classification", cost = taxaErro)
        
        # Teste
        x = extrairMatrizSemColunaDeClasse(conjuntoTestes, indiceColunaClasse)
        y = extrairColunaDeClasse(conjuntoTestes, indiceColunaClasse)
        resultadoTestes = table(predict(modelo, x), y)
        
        # Incremento dos valores de medidas
        subtotal = length(conjuntoTestes[,1])
        tpParcial = resultadoTestes[1, 1]
        fpParcial = resultadoTestes[1, 2]
        fnParcial = resultadoTestes[2, 1]
        tnParcial = resultadoTestes[2, 2]
        pParcial = tpParcial + fnParcial
        nParcial = tnParcial + fpParcial
        
        tp = tp + tpParcial
        fp = fp + fpParcial
        fn = fn + fnParcial
        tn = tn + tnParcial
        total = total + subtotal
        
        taxaFalsaAceitacao = taxaFalsaAceitacao + (fpParcial / subtotal)
        taxaFalsaRejeicao = taxaFalsaRejeicao + (fnParcial / subtotal)
        sensibilidade = sensibilidade + (tpParcial / pParcial)
        especificidade = especificidade + (tnParcial / nParcial)
        erroMedio = erroMedio + ((fpParcial + fnParcial) / subtotal)
    }
    
    # Cálculos das medidas
    p = tp + fn
    n = tn + fp
    taxaFalsaAceitacao = taxaFalsaAceitacao / k
    taxaFalsaRejeicao = taxaFalsaRejeicao / k
    sensibilidade = sensibilidade / k
    especificidade = especificidade / k
    erroMedio = erroMedio / k
    
    print(paste("Positivos:", p))
    print(paste("Negativos:", n))
    print(paste("Verdadeiros Positivos:", tp))
    print(paste("Falsos Positivos:", fp))
    print(paste("Falsos Negativos:", fn))
    print(paste("Verdadeiros Negativos:", tn))
    print(paste("Erro médio:", erroMedio))
    print(paste("Taxa de Falsa Aceitação:", taxaFalsaAceitacao))
    print(paste("Taxa de Falsa Rejeição:", taxaFalsaRejeicao))
    print(paste("Sensibilidade:", sensibilidade))
    print(paste("Especificidade:", especificidade))
}

# Dados para testes

# Matriz de testes:
# 1, 2, 3, 4
# 5, 6, 7, 8
# 9, 10, 11, 12
# 13, 14, 15, 16
# 17, 18, 19, 20
# 21, 22, 23, 24

dadosTestados = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24), 6, 4, byrow = TRUE)
numeroGrupos = 5

print(">> Início dos testes 'Dividir em K Grupos' <<")

gruposObtidos = dividirEmKGrupos(dados = dadosTestados, k = numeroGrupos)
indiceGrupo = 1

for (k.grupo in gruposObtidos) {
    print(paste("Grupo", indiceGrupo))
    print(k.grupo@grupo)
    indiceGrupo = indiceGrupo + 1
}

print(">> Fim dos testes 'Dividir em K Grupos' <<")

print(">> Início dos testes 'Avaliar SVM' <<")

# Carga da biblioteca que contém a implementação de SVM
library("e1071") 

# Carga do arquivo de dados
f = file("D:\\USP\\PPgSI - Mestrado\\2010-2_Sem\\SIN 5007-1 - Reconhecimento de Padrões\\Atividades\\Atividade 4\\Planilhas\\census_income_com_classe.csv")
censusIncomeDataSet = read.table(f, sep = ",", colClasses = c("integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer"), header = TRUE)

censusIncomeDataSetClass1 = subset(censusIncomeDataSet, Income == 1)
censusIncomeDataSetClass2 = subset(censusIncomeDataSet, Income == 2)
dataSetTeste = rbind(censusIncomeDataSetClass1[1:500,], censusIncomeDataSetClass2[1:500,])

# Teste da função de avaliação
print("Dataset: Census Income")

# Testes com diferentes valores de C (melhor C = 10)
taxasErro = c(10, 500, 2^70)

for (t in taxasErro) {
    avaliarSVM(dados = dataSetTeste, k = 7, indiceColunaClasse = 15, kernel = "radial", taxaErro = t)	
}

taxasErro = c(10, 50, 100) # melhor C = 10

for (t in taxasErro) {
    avaliarSVM(dados = dataSetTeste, k = 7, indiceColunaClasse = 15, kernel = "radial", taxaErro = t)   
}

# Testes com diferentes funções Kernel (melhor Kernel = linear)
funcoesKernel = c("linear", "radial", "polynomial", "sigmoid")

for (k in funcoesKernel) {
    avaliarSVM(dados = dataSetTeste, k = 7, indiceColunaClasse = 15, kernel = k, taxaErro = 10)   
}

# Testes com diferentes valores de K (melhor K = 10; não rodou para 999)
valoresK = c(2, 10, 100, 999)

for (v in valoresK) {
    avaliarSVM(dados = dataSetTeste, k = v, indiceColunaClasse = 15, kernel = "linear", taxaErro = 10)   
}

valoresK = c(10, 20, 50, 300, 500) # melhor K = 10 

for (v in valoresK) {
    avaliarSVM(dados = dataSetTeste, k = v, indiceColunaClasse = 15, kernel = "linear", taxaErro = 10)   
}

print(">> Fim dos testes 'Avaliar SVM' <<")

### Atividade 6 ###

# Carga do arquivo de dados
f = file("D:\\USP\\PPgSI - Mestrado\\2010-2_Sem\\SIN 5007-1 - Reconhecimento de Padrões\\Atividades\\Atividade 6\\census_income_colunas-5_6-selecionadas.csv")
dataSetTeste = read.table(f, sep = ",", colClasses = c("integer", "integer", "integer"), header = TRUE)

# Teste da função de avaliação
print("Dataset: Census Income")

# Execução
avaliarSVM(dados = dataSetTeste, k = 10, indiceColunaClasse = 3, kernel = "linear", taxaErro = 10)