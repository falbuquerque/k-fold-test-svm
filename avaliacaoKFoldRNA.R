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
# Avalia o classificador baseado em RNA utilizando a técnica de cross-validation
#
avaliarRNA = function(dados, k, indiceColunaClasse, numeroCamadasInternas = 1, taxaAprendizado = 0.1) {
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
    print(paste("Número de camadas internas:", numeroCamadasInternas))
    print(paste("Taxa de aprendizado:", taxaAprendizado))
    
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
        modelo = nnet(x, y, size = numeroCamadasInternas, decay = taxaAprendizado)
        
        # Teste
        x = extrairMatrizSemColunaDeClasse(conjuntoTestes, indiceColunaClasse)
        y = extrairColunaDeClasse(conjuntoTestes, indiceColunaClasse)
      
        pred = predict(modelo, x)
        
        # Somente uma coluna é retornada, com os valores aproximados para as classes
        for (x in 1:length(pred[,1])) {
            pred[x,1] = round(pred[x,1], 0)
        }

        # Incremento dos valores de medidas
        resultadoTestes = table(y, pred)

        subtotal = length(conjuntoTestes[,1])
        tpParcial = resultadoTestes[1, 1]
        fnParcial = resultadoTestes[2, 1]
        
        if (length(resultadoTestes[1,]) > 1) {
            fpParcial = resultadoTestes[1, 2]
            tnParcial = resultadoTestes[2, 2]
        } else {
            fpParcial = 0
            tnParcial = 0
        }

        total = total + subtotal
        
        pParcial = tpParcial + fnParcial
        nParcial = tnParcial + fpParcial

        tp = tp + tpParcial
        tn = tn + tnParcial
        fp = fp + fpParcial
        fn = fn + fnParcial
        
        taxaFalsaAceitacao = taxaFalsaAceitacao + (fpParcial / subtotal)
        taxaFalsaRejeicao = taxaFalsaRejeicao + (fnParcial / subtotal)
        
        if (pParcial > 0) {
            sensibilidade = sensibilidade + (tpParcial / pParcial)    
        }
        
        if (nParcial > 0) {
            especificidade = especificidade + (tnParcial / nParcial)	
        }
        
        erroMedio = erroMedio + ((fpParcial + fnParcial) / subtotal)
    }
    
    # Cálculos das medidas
    p = tp + fn
    n = tn + fp
    
    taxaFalsaAceitacao = taxaFalsaAceitacao / k
    taxaFalsaRejeicao = taxaFalsaRejeicao / k
    sensibilidade = sensibilidade / k
    
    if (especificidade > 0) {
        especificidade = especificidade / k	
    }
    
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

print(">> Início dos testes 'Avaliar Redes Neurais' <<")

# Carga da biblioteca que contém a implementação de RNA
library("nnet") 

# Carga do arquivo de dados
f = file("D:\\USP\\PPgSI - Mestrado\\2010-2_Sem\\SIN 5007-1 - Reconhecimento de Padrões\\Atividades\\Atividade 6\\census_income_colunas-5_6-selecionadas.csv")
dataSetTeste = read.table(f, sep = ",", colClasses = c("integer", "integer", "integer"), header = TRUE)

# Teste da função de avaliação
print("Dataset: Census Income")

# Testes com diferentes números de camadas (todos retornaram valores iguais)
numerosCamadas = c(1, 10, 30)

for (n in numerosCamadas) {
    avaliarRNA(dados = dataSetTeste, k = 10, indiceColunaClasse = 3, numeroCamadasInternas = n)  
}

# Testes com diferentes valores de taxa de aprendizado (todos retornaram valores iguais)
taxas = c(1e-10, 1e-2, 1, 1e2, 1e10)

for (t in taxas) {
    avaliarRNA(dados = dataSetTeste, k = 10, indiceColunaClasse = 3, numeroCamadasInternas = 10, taxaAprendizado = t)
}

#Testes com todas as características (todos retornaram valores iguais)
f = file("D:\\USP\\PPgSI - Mestrado\\2010-2_Sem\\SIN 5007-1 - Reconhecimento de Padrões\\Atividades\\Atividade 6\\census_income_com_classe_amostra.csv")
dataSetTeste = read.table(f, sep = ";", colClasses = c("integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer"), header = TRUE)

avaliarRNA(dados = dataSetTeste, k = 10, indiceColunaClasse = 15, numeroCamadasInternas = 10, taxaAprendizado = 0.01)