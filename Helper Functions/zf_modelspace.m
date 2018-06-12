function models = zf_modelspace

mat         = zeros(10,10);
dia         = mat;
for m = 1:length(mat), dia(m,m) = 1;  end  

Hom         = mat;
for m = 2:length(mat)
    Hom(m, m-1) = 1;
    Hom(m-1, m) = 1;
end

Nei         = mat;
for m = 3:length(mat)
    Nei(m, m-2) = 1;
    Nei(m-2, m) = 1;
end

Hublbls = {'Tec', 'Cbl', 'RHb', 'CHb', 'RSC'};

m = 0;
clear M; 

for n = 1:2
for l = 1:2
if n == 1, nei = mat; nstr = ''; else nei = Nei; nstr = '_nei'; end
if l == 1, hom = mat; lstr = ''; else hom = Hom; lstr = '_hom'; end
    
for h = 0:5
    m = m+1;

    M{m}.A{1}   = mat + triu(nei) + hom;
    M{m}.A{2}   = mat + tril(nei) + hom;
    M{m}.A{3}   = dia;
    
    if h ~= 0  
        M{m}.name = ['Hub_' Hublbls{h} nstr lstr];
        for a = 1:length(M{1}.A{1})
           if (1+2*(h-1)) < a, 
                M{m}.A{1}(1+2*(h-1),a) = 1;
                M{m}.A{2}(a,1+2*(h-1)) = 1;
           elseif (1+2*(h-1)) > a, 
                M{m}.A{2}(1+2*(h-1),a) = 1;
                M{m}.A{1}(a,1+2*(h-1)) = 1;
           end

           if (2+2*(h-1)) < a, 
                M{m}.A{1}(2+2*(h-1),a) = 1;
                M{m}.A{2}(a,2+2*(h-1)) = 1;
           elseif (2+2*(h-1)) > a, 
                M{m}.A{2}(2+2*(h-1),a) = 1;
                M{m}.A{1}(a,2+2*(h-1)) = 1;
           end
        end
    else
        M{m}.name = ['Hub_None' nstr lstr];
    end
end

end
end

lm      = length(M) + 1;
M{lm}.name   = 'full';
ful          = ones(10,10) - dia;
M{lm}.A{1}   = ful;
M{lm}.A{2}   = ful;
M{lm}.A{3}   = mat + dia;

models = M;
