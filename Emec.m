% Nom du fichier d'output a analyser (modifiez selon vos besoins)
filename = 'nsteps=10000.out'; 

tfin   = 1e-6;      % temps total maximal de simulation [s]
nsteps   = 10000;
dt = tfin/nsteps;
% Chargement des donnees
data = load(filename);

t = 0;
ListeTemps = [];
while tfin > t
    ListeTemps = [ListeTemps; t];
    t = t + dt;
end

% Extraction des quantites d'interet
Emeca = data(:,6);
Emec0 = data(1,6);

%Calcul de la pente de Emec par rapport à t
PenteEmec = [];
for c = 1:(size(Emeca)-1)
    PenteEmec = [PenteEmec; (Emeca(c)-Emec0)];
end

size(PenteEmec)
size(ListeTemps)

%Génération du graphe
lw=2; fs=16;
figure
plot(ListeTemps, PenteEmec, 'k+-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('Temps [s]')
ylabel('\Delta E_{mec} [J]')
grid on
