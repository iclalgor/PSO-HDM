function y = hiperkureUzerindeOtelemenew(M,r,dphi,w,yon)
% n-boyutlu uzayda M merkezli r yarıçaplı hiperküre üzerinde yer alan P
% noktasını hiperküreyi ilk eksen ektarında dphi açısı kadar döndürülünye
% yeni noktanın koordinatlarını veren fonksiyon
% M - Hiperkürenin merkez koordinatları
% r - Hiperkürenin yarıçapı
% dphi - Döndürme açısı
% w - Hiperküre üzerindeki bir nokta
% yon - Eğer yon -1 ise saat yönünün tersi
% Son Güncelleme : 16.11.2012
w=w-M;
n = size(M,2);%numel(M.Position);%10;%size(M,2); % Uzayın boyutu

% Eğer n = 1 ise, hiper küre tanımı geçerli değil, bu yüzden hata versin
% if n < 2
%     error('Hiper küre en az iki boyutlu bir uzayda tanımlıdır.');
% end

% w noktası hiper küre üzerinde olduğu için denklemi sağlar. O halde açı
% değerlerini hesaplatabiliriz.
phi=zeros(1,n-1);
for i=1:n-2
    toplam=0;
    for j=i+1:n
        toplam=toplam+w(j)^2;
    end
    phi(i)=acot(w(i)/sqrt(toplam));
end
 phi(n-1)=2*acot((sqrt(w(n-1)^2+w(n)^2)+w(n-1))/w(n));

%% w noktasını ilk eksen yönünde dphi kadar döndürürürsek yeni koordinatlar
% ne olur?
phi(1)=phi(1)-yon*(phi(1)+dphi);

w(1)=r*cos(phi(1));
c=1;
for i=2:n-1
    c=1;
    for k=1:i-1
        c = c*sin(phi(k));
    end
    w(i) = r*c*cos(phi(i));
end
w(n) =r*c*sin(phi(n-1));

%Nokta döndürme işlemi tamamlandıktan sonra tekrar öteleniyor.
w=w+M;
y=w;
end %function

