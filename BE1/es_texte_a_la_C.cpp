#ifdef IO_TEXT_C

#include <cstdio> //fopen, fclose, fprintf, fscanf
#include <cstdlib> //exit

void creeFichierEntiers(const char * nomFichier, int nb)
//preconditions : nomFichier chaine de caracteres designant le nom du fichier a creer
//postcondition : le fichier nomFichier contient nb entiers separes par des espaces
{
  FILE *fd=std::fopen(nomFichier,"w");
  if (fd==NULL)
    {std::printf("Impossible d'ouvrir le fichier %s en ecriture \n",nomFichier);exit(1);}

  for(int i=0;i<nb;i++)
    {
      int temp;
      std::printf("Entier suivant : ");
      std::scanf("%d",&temp); //Ou utilisez la version robuste de saisie d'un int
      std::fprintf(fd,"%d ",temp);//Remarquez que l'on separe les int par des espaces
    }
  std::fclose(fd);
}

void litFichierEntiers(const char * nomFichier)
//preconditions : nomFichier chaine de caracteres designant le nom d'un fichier
//                contenant des entiers separes par des caracteres d'espacement
//postcondition : affichage du contenu du fichier sur la sortie standard
{
  FILE *fd=std::fopen(nomFichier,"r");
  if (fd==NULL)
    {std::printf("Impossible d'ouvrir le fichier %s en lecture \n",nomFichier);exit(1);}

  int temp;
  while(std::fscanf(fd,"%d",&temp)==1)
    std::printf("Entier suivant : %d \n",temp);
  std::fclose(fd);
}

int main()
{
  char fich[100]; 
  int n;
  std::printf("Entrez le nom du fichier d'entiers a creer\n");
  std::scanf("%99s",fich);
  std::printf("Combien d'entiers dans ce fichier ?\n");
  std::scanf("%d",&n);
  creeFichierEntiers(fich,n);
  std::printf("Affichage des entiers contenu dans le fichier\n");
  litFichierEntiers(fich);
  return 0;
}

#endif
