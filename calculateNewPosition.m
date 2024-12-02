function newPos = calculateNewPosition(rho, alpha, h, nVar)
    newPos = zeros(1, nVar);  % Yeni pozisyonu başlat
    newPos(1) = rho * cos(alpha(1) + h);  % İlk boyut

    % Diğer boyutlar
    for ii = 2:nVar-1
        c = 1;
        for k = 1:ii-1
            c = c * sin(alpha(k));  % Çarpım sinüs terimlerini al
        end
        newPos(ii) = rho * c * cos(alpha(ii) + h);  % Diğer boyutların pozisyonunu hesapla
    end
    % Son boyut için c tüm önceki alpha'ların sinüslerinin çarpımı olmalı
    c = 1;
    for k = 1:nVar-1
        c = c * sin(alpha(k));
    end
    newPos(nVar) = rho * c * sin(alpha(nVar-1) + h);  % Son boyut
end
