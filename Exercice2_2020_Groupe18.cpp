#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
                          // Fichier .tpp car inclut fonctions template

using namespace std; // ouvrir un namespace avec la librerie c++ de base

/* definir a fonction template pour calculer le produit interne
   entre deux valarray
   inputs:
     array1: (valarray<T>)(N) vecteur de taille N
     array2: (valarray<T>)(N) vecteur de taille N
   outputs:
     produitInterne: (T) produit entre les deux vecteurs
*/
template<typename T> T scalarProduct(valarray<T> const& array1,\
valarray<T> const& array2){
  // compute and return the norm2 of a valarray
  return (array1*array2).sum();
} 

/* definir a fonction template pour calculer la norm2 d'un valarray
   inputs:
     array: (valarray<T>)(N) vecteur de taille N
   outputs:
     norm2: (T) norm2 du vecteur
*/
template<typename T> T norm2(valarray<T> const& array){
  // compute and return the norm2 of a valarray
  return sqrt((array*array).sum());
} 

/* definir a fonction template pour calculer le produit vecteur
   entre 2 valarray de dimension 3
   inputs:
     array1, array2: (valarray<T>)(N) vecteurs de taille N
   outputs:
     produitVecteur: (T) produit vectoriel array1 x aray2 
*/
template<typename T> T vectorProduct(valarray<T> const& array1,\
valarray<T> const& array2){
  // initialiser le nouveau valarray
  valarray<T> array3=valarray<T>(3);
  // calculer le produit vecteur
  array3[0] = array1[1]*array2[2] - array1[2]*array2[1]; // premier composante
  array3[1] = array1[2]*array2[0] - array1[0]*array2[2]; // deuxieme composante
  array3[2] = array1[0]*array2[1] - array1[1]*array2[0]; // troisieme composante
  // retourner le array3
  return array3;
} 


/* La class Engine est le moteur principale de ce code. Il contient 
   les methodes de base pour lire / initialiser les inputs, 
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{

private:
  // definition des constantes
  const double pi=3.1415926535897932384626433832795028841971e0;
  // definition des variables
  double tfin=0.e0;      // Temps final
  unsigned int nsteps=1; // Nombre de pas de temps
  double charge=0.e0; 	 // charge de la particule
  double mass=1.e0;      // mass de la particule
  double B0=0.e0;	 // intensite du champ magnetique
  double l_k=0.e0;       // longueur caracteristique du champ magnetique 
  double E0=0.e0;	 // intensite du champ electrique
  valarray<double> x0=valarray<double>(0.e0,2); // vecteur contenant la position initiale du ballon en
			 		        // utilisant la sequence index-valeur: 0-x, 1-z
  valarray<double> v0=valarray<double>(0.e0,2); // vecteur contenant la vitesse initiale du ballon en
			  		        // utilisant la sequence index-valeur: 0-vx, 1-vz
  unsigned int sampling=1; // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  // TODO calculer l'energie mecanique et le moment magnetique
  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
       write: (bool) ecriture de tous les sampling si faux
  */  
  void printOut(bool write)
  {
    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      double mechanicalEnergy = 0.5*mass*(v[0]*v[0] + v[1]*v[1]) - charge*E0*x[1]; // TODO calculer l'energie mecanique 
      double magneticMoment   = mass*(v[0]*v[0] + v[1]*v[1])/(2*B0); // TODO calculer le moment magnetique 
      *outputFile << t << " " << x[0] << " " << x[1] << " " \
      << v[0] << " " << v[1] << " " << mechanicalEnergy << " " \
      << magneticMoment << endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

  // Iteration temporelle, a definir au niveau des classes filles
  virtual void step()=0;

protected:

  // donnes internes
  double t,dt;        // Temps courant pas de temps
  double chargeOmass; // charge / mass
  valarray<double> x=valarray<double>(2); // Position actuelle de la particule 
  valarray<double> v=valarray<double>(2); // Vitesse actuelle de la particule

   /* Cette methode calcule le champ magnetique en y
     outputs:
       By: (double) champ magnetique en y
  */ 
  double By() const
  {
     return B0*(1.e0+l_k*x[0]);
  }

  /* Cette methode calcule le champ electrique en z
     outputs:
       Ez: (double) champ electrique en z
  */ 
  valarray<double> Ez() const
  {
    valarray<double> E=valarray<double>({0.e0,E0}); 
    return E;
  }

  // TODO
  /* Cette methode calcule l'acceleration
     outputs:
       a: (valarray<double>)(2) acceleration dans les directions (x,z)
  */
  valarray<double> acceleration(valarray<double>& a) const
  {
    // compute the acceleration
    a[0] = -(charge/mass)*By()*v[1]; // TODO
    a[1] = (charge/mass)*(By()*v[0] + Ez()[1]); // TODO
    return a;
  }

public:

  /* Constructeur de la classe Engine
     inputs:
       configFile: (ConfigFile) handler du fichier d'input
  */
  Engine(ConfigFile configFile)
  {
    // variable locale
    
    // Stockage des parametres de simulation dans les attributs de la classe
    tfin     = configFile.get<double>("tfin");		 // lire la temps totale de simulation
    nsteps   = configFile.get<unsigned int>("nsteps");   // lire la nombre de pas de temps
    mass     = configFile.get<double>("mass");		 // lire la mass de la particule
    charge   = configFile.get<double>("charge");	 // lire la charge de la particule
    B0       = configFile.get<double>("B0");		 // lire l intensite du champ magnetique
    l_k      = configFile.get<double>("l_k");		 // lire la longueur du gradient de B
    E0       = configFile.get<double>("E0");		 // lire l intensite du champ electrique
    x0[0]    = configFile.get<double>("x0");		 // lire composante x position initiale
    x0[1]    = configFile.get<double>("z0");		 // lire composante z position initiale
    v0[0]    = configFile.get<double>("vx0");		 // lire composante x vitesse initiale
    v0[1]    = configFile.get<double>("vz0");		 // lire composante z vitesse initiale
    sampling = configFile.get<unsigned int>("sampling"); // lire le parametre de sampling

    dt = tfin / nsteps;          // calculer le time step
    chargeOmass = charge / mass; // coefficient charge/mass

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str()); 
    outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales
  };

  // Destructeur virtuel
  virtual ~Engine()
  {
    outputFile->close();
    delete outputFile;
  };

  // Simulation complete
  void run()
  {
    t = 0.e0; // initialiser le temps
    x = x0;   // initialiser la position
    v = v0;   // initialiser la vitesse
    last = 0; // initialise le parametre d'ecriture
    printOut(true); // ecrire premier pas de temps
    for(unsigned int i(0); i<nsteps; ++i) // boucle sur tout pas de temps
    {
      step();  // faire la mise a jour de la simulation 
      printOut(false); // ecrire pas de temps actuel
    }
    printOut(true); // ecrire dernier pas de temps
  };

};

// Extension de la class Engine implementant l'integrateur d'Euler
class EngineEulerImplicite: public Engine
{
private:
  unsigned int maxit=1000; // nombre maximale d iterations
  double tol=1.e12;        // tolerance methode iterative
public:
  
  // construire la class Engine
  EngineEulerImplicite(ConfigFile configFile): Engine(configFile){

    tol = configFile.get<double>("tol"); // lire la tolerance pour la method iterative
    maxit = configFile.get<unsigned int>("maxit"); // lire le nombre d iterations maximale 

  }

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le schema: Euler-Implicit en utilisant la methode de point fixe
     pour resoude le problem implicite
  */
  void step()
  {
    unsigned int iteration=0;
    //double error=999e0;
    valarray<double> a=valarray<double>(0.e0,2);
    acceleration(a);
    valarray<double> xold=valarray<double>(x);
    valarray<double> vold=valarray<double>(v);
    double norme(0);
    valarray<double> ErreurTemp(0.0, 4);
    
    do{
		x = xold + v*dt;
		v = vold + a*dt;
		acceleration(a);
		
		ErreurTemp = {x[0] - xold[0] - v[0]*dt, x[1] - xold[1] - v[1]*dt, v[0] - vold[0] - a[0]*dt, v[1] - vold[1] - a[1]*dt}; // Concaténation des vecteurs position et vitesse pour calculer la norme
		norme = norm2(ErreurTemp);
		++iteration;
	} while((norme > tol) and (iteration<=maxit));
    
    t += dt; // mis a jour du temps
  }
};

// Extension de la class Engine implementant l'integrateur d'Euler
class EngineEuler: public Engine
{
public:
  
  // construire la class Engine
  EngineEuler(ConfigFile configFile): Engine(configFile) {}

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le schema: Euler
  */
  void step()
  {
    valarray<double> a=valarray<double>(0.e0,2);
    acceleration(a);
    
    x += v*dt;
    v += a*dt;
    t += dt; // mis a jour du temps
  }
};

// Extension de la class Engine implementant l'integrateur d'Euler-Cromer
class EngineEulerCromer: public Engine
{
public:

  // construire la class Engine
  EngineEulerCromer(ConfigFile configFile): Engine(configFile) {}

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le scheme: Euler-Cromer
  */
  void step()
  {
     valarray<double> a=valarray<double>(0.e0,2);
     acceleration(a);
     
	 x+=v*dt;
	 v[0]+=a[0]*dt;
	 acceleration(a);
	 v[1]+=a[1]*dt;
	 
     t += dt; // mis a jour du temps
  }
};

// Extension de la class Engine implementant l'integrateur Runge-Kutta 2
class EngineRungeKutta2: public Engine
{
public:

  // construire la class Engine
  EngineRungeKutta2(ConfigFile configFile): Engine(configFile) {}

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le scheme: Runge-kutta 2
  */
  void step()
  {
    valarray<double> a=valarray<double>(0.e0,2);
    //valarray<double> xold=valarray<double>(x);
    //valarray<double> vold=valarray<double>(v);
    acceleration(a); // definir l'acceleration
    
    x += 0.5*v*dt; // t_(i+1/2)
    v += 0.5*a*dt; // t_(i+1/2)
    acceleration(a); // Mise à jour (t_(i+1/2))
    v += 0.5*a*dt; // neuf
    x += 0.5*v*dt; // neuf
    //x = xold + v*dt;
    //v = vold + a*dt;

    t += dt; // mis a jour du temps   
  }
};


// Extension de la class Engine implementant l'integrateur Runge-Kutta 2
class EngineBorisBuneman: public Engine
{

protected:

  // TODO
  // Calcul la rotation des vitesses pour la method de Boris Buneman
  void rotationVitessesBorisBuneman()
  {
      double w_c = chargeOmass*By();
	  double b = w_c*dt;
	  
	  double vx = ((1-pow(b,2)/4)*v[0]-b*v[1]) / (1+pow(b,2)/4);
	  double vz = ((1-pow(b,2)/4)*v[1]+ b*v[0]) / (1+pow(b,2)/4);
	  
	  v[0] = vx;
	  v[1] = vz;
  }

public:

  // construire la class Engine
  EngineBorisBuneman(ConfigFile configFile): Engine(configFile) {}

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le scheme: Boris-Buneman. 
  */
  void step()
  {
	valarray<double> a=valarray<double>(0.e0,2);
    acceleration(a);
    
	x += v*dt/2; 
    v += chargeOmass*Ez()*dt/2; 
    
    rotationVitessesBorisBuneman();
    
    v += chargeOmass*Ez()*dt/2; 
    x += v*dt/2;

    t += dt; // mis a jour du temps
  }
};

// programme
int main(int argc, char* argv[])
{
  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Schema numerique ("Euler"/"E", "EulerCromer"/"EC" ou "RungeKutta2"/"RK2")
  string schema(configFile.get<string>("schema"));

  Engine* engine; // definer la class pour la simulation
  // choisir quel schema utiliser
  if(schema == "Euler" || schema == "E")
  {
    // initialiser une simulation avec schema Euler
    engine = new EngineEuler(configFile);
  }
  else if(schema == "EulerImplicite" || schema == "EI")
  {
    // initialiser une simulation avec schema Euler Implicite
    engine = new EngineEulerImplicite(configFile);
  }
  else if(schema == "EulerCromer" || schema == "EC")
  {
    // initialiser une simulation avec schema Euler-Cromer
    engine = new EngineEulerCromer(configFile);
  }
  else if(schema == "RungeKutta2" || schema == "RK2")
  {
    // initialiser une simulation avec schema runge-kutta 2
    engine = new EngineRungeKutta2(configFile);
  }
  else if(schema == "BorisBuneman" || schema == "BB")
  {
    // initialiser une simulation avec schema Boris Buneman
    engine = new EngineBorisBuneman(configFile);
  }
  else
  {
    cerr << "Schema inconnu" << endl;
    return -1;
  }

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation 
  cout << "Fin de la simulation." << endl;
  return 0;
}
