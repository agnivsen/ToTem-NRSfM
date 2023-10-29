addpath('D:\code\SEDUMI\sedumi-master\');
addpath('D:\code\SfT-Experiments_v2\matlabscripts\Local-Iso-NRSfM_v1p1\gloptipoly3\');


P = zeros(12,12);

P(:,1) = [1 2 5 4 8 7 6 3 5 7 4 3].';


C{1}.c = P;

C{1}.t = 'min';

gloptipoly(C);