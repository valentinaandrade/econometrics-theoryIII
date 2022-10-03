function [Tendencia,Ciclo] = HP(Y,smoothing)
% Referencia: Hodrick, R. J., and E. C. Prescott (1997). "Postwar U.S. Business Cycles:
% An Empirical Investigation." Journal of Money, Credit, and Banking.Vol. 29, pp. 1-16.
% 
%% Detalles de la Función HP:
%
% HP: Filto Hodrick-Prescott estima los componentes de tendencia y ciclo de
% la o las series de interes
%
% Sintaxis:
%	[Tendencia,Ciclo] = HP(Y,smoothing)
%	 HP(...)
%
% Descripción:
%   Separe una o más series de tiempo en componentes de tendencia y ciclo 
%   con el filtro de Hodrick-Prescott. Si no se especifica argumentos de salida,
%   la función HP muestra un gráfico de la serie y la tendencia (con los ciclos 
%   eliminados). El gráfico se puede utilizar para ayudar a seleccionar un parámetro
%   de suavizamiento.
%
% Input (Series):
%   Y - Datos de series de tiempo. puede ser un vector o una matriz. Si Y es 
%       un vector, representa una sola serie. Si es una matriz n x m, representa 
%       las n observaciones de las m series, y se supone que las observaciones 
%       en cualquier fila ocurren al mismo tiempo. Además, se supone que la última 
%       observación de cualquier serie es la más reciente.
%
% Input Adicional (Factor de Suavizamiento / smoothing):
%   smoothing: un escalar que se aplicará a todas las series o un vector de
%       longitud n con valores que se aplicarán a las series correspondientes.
%       El valor predeterminado es 1600, que se sugiere para datos trimestrales. 
%       Si el smoothing=0, no se produce ningún suavizado. A medida que aumenta el 
%       parámetro de suavizado, la serie suavizada se aproxima a una línea recta. 
%       Si el suavizado es infinito (Inf), la serie no tiene tendencia.
%
% Output:
%	Tendencia - Componente de tendencia de Y, tiene el mismo tamaño que Y.
%	Ciclo     - Componenete ciclico de Y, tiene el mismo tamaño que Y.
%
% Observaciones: 
%       [1] Y = Tendencia + Ciclo
%       [2] Valores sugeridos para el parámetro smoothing
%               - Frecuencia Anual = 100
%               - Frecuencia Trimestral = 1600
%               - Frecuencia Mensual = 14400
%       [3] El filtro de Hodrick-Prescott puede producir sesgo del punto final 
%           en datos de muy alta frecuencia.

%% Función HP

% Comprobamos los argumentos de entrada

if nargin < 1 || isempty(Y)
	error('FaltanDatosEntrada')
end

if ~isscalar(Y) && isvector(Y) && isa(Y,'double')
	Y = Y(:);
	[numObs,numSeries] = size(Y);
elseif ndims(Y) == 2 && min(size(Y)) > 1 && isa(Y,'double')
	[numObs,numSeries] = size(Y);
else  
	error('EtradaNoValidaArg1')
end

if any(any(~isfinite(Y)))
	error('DatosEntradaNoValidos') 
end

if numObs < 3 % Tratar las muestras con < 3 observaciones solo como datos de tendencia
	warning('DatosInsuficientes')
	Tendencia = Y;
	Ciclo = zeros(numObs,numSeries);
	return  
end

if nargin < 2 || isempty(smoothing)
	warning('TrimPredetSmoothing')
	smoothing = 1600;  
end

if ~isvector(smoothing) || ~isa(smoothing,'double')
	error('EntradaNoValidaArg2')
    
else
	if ~any(numel(smoothing) == [1,numSeries])
		error('DimesnionSmoothingIncosis')
	end
end

if any(isnan(smoothing))  
	error('SmoothingInvalido')
end

if any(smoothing < 0)
	warning('SmoothingNegativo')
	smoothing = abs(smoothing);
end

% Ejecute el filtro con factor smoothing escalar o vectorial:

if (numel(smoothing) == 1) || (max(smoothing) == min(smoothing)) % Parámetro smoothing
    if numel(smoothing) > 1
		smoothing = smoothing(1);
    end
	if isinf(smoothing)	% Usa MCO para quitar la tendencia        
		Tendencia = Y-detrend(Y);    
    else    
		if numObs == 3 % Caso especial con 3 muestras 
			A = eye(numObs,numObs) + ...
				smoothing*[ 1 -2 1; -2 4 -2; 1 -2 1 ];    
        else % Caso generial con > 3 muestras     
			e = repmat([smoothing,-4*smoothing,(1+6*smoothing),...
				        -4*smoothing,smoothing],numObs,1);
			A = spdiags(e,-2:2,numObs,numObs);
			A(1,1) = 1+smoothing;
			A(1,2) = -2*smoothing;
			A(2,1) = -2*smoothing;
			A(2,2) = 1+5*smoothing;
			A(numObs-1,numObs-1) = 1+5*smoothing;
			A(numObs-1,numObs) = -2*smoothing;
			A(numObs,numObs-1) = -2*smoothing;
			A(numObs,numObs) = 1+smoothing;    
        end
		Tendencia = A\Y;  
	end
    
else % Vector smoothing
    
	Tendencia = zeros(numObs,numSeries);    
	if numObs == 3 % Caso especial con 3 muestras
		for i = 1:numSeries            
            if isinf(smoothing(i)) % Use OLS detrending
				Tendencia(:,i) = Y(:,i)-detrend(Y(:,i));  
            else 
				A = eye(numObs,numObs) + ...
					smoothing(i)*[ 1 -2 1; -2 4 -2; 1 -2 1 ];
				Tendencia(:,i) = A\Y(:,i);
            end
		end
        
	else % General case with > 3 samples
        
        for i = 1:numSeries  
            if isinf(smoothing(i)) % Usa MCO para quitar la tendencia 
				Tendencia(:,i) = Y(:,i)-detrend(Y(:,i));           
            else       
				e = repmat([smoothing(i),-4*smoothing(i),(1+6*smoothing(i)), ...
					        -4*smoothing(i),smoothing(i)],numObs,1);
				A = spdiags(e,-2:2,numObs,numObs);
				A(1,1) = 1+smoothing(i);
				A(1,2) = -2*smoothing(i);
				A(2,1) = -2*smoothing(i);
				A(2,2) = 1+5*smoothing(i);
				A(numObs-1,numObs-1) = 1+5*smoothing(i);
				A(numObs-1,numObs) = -2*smoothing(i);
				A(numObs,numObs-1) = -2*smoothing(i);
				A(numObs,numObs) = 1+smoothing(i);
				Tendencia(:,i) = A\Y(:,i);   
            end
        end
    end
end

% Si no hay argumentos de salida, grafica por defecto los siguientes resultados:

if nargout == 0    
	figure(gcf);
	plot(Y,'b');
	hold on
	plot(Tendencia,'r');
	title('\bfFiltro Hodrick-Prescott');
    
	if numSeries == 1
		legend('Serie','Tendencia Suavizada');
    end
	hold off;   
elseif nargout > 1   
	Ciclo = Y-Tendencia;  
end