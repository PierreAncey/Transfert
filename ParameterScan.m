% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = 'posixfs/Exercice2_2020/'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice2_2020_Groupe18'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'posixfs/Exercice2_2020/configuration.in'; % Nom du fichier d'entree de base MODIFIER SELON VOS BESOINS

nsimul = 5; % Nombre de simulations a faire

nsteps = round(logspace(2,4,nsimul)); % Nombre d'iterations entier de 10^2 a 10^4  MODIFIER SELON VOS BESOINS

paramstr = 'nsteps'; % Nom du parametre a scanner  MODIFIER SELON VOS BESOINS
param = nsteps; % Valeurs du parametre a scanner  MODIFIER SELON VOS BESOINS

%% Simulations %% 
%%%%%%%%%%%%%%%%%
% Lance une serie de simulations (= executions du code C++)
% Normalement, on ne devrait pas avoir besoin de modifier cette partie

output = cell(1, nsimul); % Tableau de cellules contenant les noms des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i)), '.out']
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i})
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%
% Ici, on aimerait faire une etude de convergence: erreur fonction de dt, sur diagramme log-log.
% A MODIFIER ET COMPLETER SELON VOS BESOINS

tfin   = 1e-8;      % temps total maximal de simulation [s]
mass   = 1.6726e-27;  % mass de la particule test [kg]
charge = 1.6022e-19;  % charge de la particule test [C] 
B0     = 5e0; 	     % intensite du champ magnetique [T]
l_k    = 0e0;	     % longueur characteristique due champ magnetique [1/m]
E0     = 0e0;        % intensite du champ electrique [V/m]
x0     = 0e0;        % position initiale en x [m]
z0     = 0e0;        % position initiale en z [m]
vx0    = 0e0;        % vitesse initiale en x [m/s]
vz0    = 4e5;        % vitesse initiale en z [m/s]
w0 = B0*charge/mass;

errorpos = zeros(1,nsimul);
errorspeed = zeros(1,nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    dt(i) = data(2,1)-data(1,1); % Delta t = t_2-t_1
    xend = data(end,2); % Extraire le x final
    zend = data(end,3); % Extraire le z final
    vxend = data(end,4); % Extraire le x final
    vzend = data(end,5); % Extraire le z final
    x_th = x(tfin, E0, B0, w0, vz0); % TODO: Entrer la vraie solution analytique a tfin
    z_th = z(tfin, E0, B0, w0, vz0); % TODO: Entrer la vraie solution analytique a tfin
    vx_th = vx(tfin, E0, B0, w0, vz0); % TODO: Entrer la vraie solution analytique a tfin
    vz_th = vz(tfin, E0, B0, w0, vz0); % TODO: Entrer la vraie solution analytique a tfin
    errorpos(i) = sqrt((xend-x_th)^2+(zend-z_th)^2); % erreur sur la position finale
    errorspeed(i) = sqrt((vxend-vx_th)^2+(vzend-vz_th)^2); % erreur sur la position finale
end
%Errorpos =max(errorpos)
%Errorspeed =max(errorspeed)

%figure
%loglog(dt, errorpos, 'k+')
%xlabel('\Delta t')
%ylabel('Erreur sur la position finale')
%grid on
errorpos
errorspeed

function x = x(t, E0, B0, w0, vz0)
    x = (E0/B0)*(sin(w0*t)/w0 - t) + (vz0/w0)*(cos(w0*t) - 1);
end
function z = z(t, E0, B0, w0, vz0)
    z = (E0/(B0*w0))*(1-cos(w0*t)) + (vz0/w0)*sin(w0*t);
end
function vx = vx(t, E0, B0, w0, vz0)
    vx = (E0/B0)*(cos(w0*t) - 1) - vz0*sin(w0*t);
end
function vz = vz(t, E0, B0, w0, vz0)
    vz = (E0/B0)*sin(w0*t) + vz0*cos(w0*t);
end


