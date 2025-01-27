#ifdef IO_TEXT_CPP

#include <iostream> 
#include <cstdio>
#include <fstream> //ifstream, ofstream
#include <cstdlib> //exit

void creeFichierEntiers(const char * nomFichier, int nb)
//preconditions : nomFichier chaine de caracteres designant le nom du fichier a creer
//postcondition : le fichier nomFichier contient nb entiers separes par des espaces
{
  std::ofstream ofs;
  ofs.open(nomFichier);
  if(ofs.bad()) 
    {std::cout<<"Impossible d'ouvrir le fichier "<<nomFichier<<" en ecriture \n"; exit(1);}

  for(int i=0;i<nb;i++)
    {
      int temp;
      std::cout << "Entier suivant : ";
      std::cin >> temp; //Ou utilisez la version robuste de saisie d'un int
      ofs << temp <<' ';//Remarquez que l'on separe les int par des espaces
    }
  ofs.close();
}

void litFichierEntiers(const char * nomFichier)
//preconditions : nomFichier chaine de caracteres designant le nom d'un fichier
//                contenant des entiers separes par des caracteres d'espacement
//postcondition : affichage du contenu du fichier sur la sortie standard
{
  std::ifstream ifs;
  ifs.open(nomFichier);
  if(ifs.bad())
    {std::cout<<"Impossible d'ouvrir le fichier "<<nomFichier<<" en lecture \n"; exit(1);}
  int temp;
  while(ifs >> temp)
    std::cout<<"Entier suivant "<<temp<< std::endl;
  ifs.close();
}

int main()
{
  char fich[100]; 
  int n;
  std::cout << "Entrez le nom du fichier d'entiers a creer" << std::endl;
  std::scanf("%99s",fich);
  std::cout << "Combien d'entiers dans ce fichier ?" << std::endl;
  std::cin >> n;
  creeFichierEntiers(fich,n);
  std::cout << "Affichage des entiers contenu dans le fichier" << std::endl;
  litFichierEntiers(fich);
  return 0;
}

#endif