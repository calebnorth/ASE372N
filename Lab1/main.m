
name = '3a';
input  = strcat(name,'.csv');
output = strcat(name,'.txt');
kmlput = strcat(name,'.kml');

[sod, lats, longs, hts, sats] = readNMEAposV3(input);

M = [sod; lats; longs; hts; sats];
N = reshape(M,[],5);
dlmwrite(output, N, 'delimiter', ',', 'precision', 12)

K = [longs; lats; hts];
L = reshape(K,[],3);
dlmwrite(kmlput, L, 'delimiter', ',', 'precision', 12)
