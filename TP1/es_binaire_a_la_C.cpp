#ifdef IO_BIN_C

#include <cstdio> //fopen, fclose, fread, fwrite, feof
#include <cstdlib> //exit

void creeFichierEntiersBin(const char * nomFichier, int nb)
//preconditions : nomFichier chaine de caracteres designant le nom du fichier a creer
//postcondition : le fichier nomFichier contient nb entiers ecrits sous forme binaire
{
  FILE *fd=std::fopen(nomFichier,"w");
  if (fd==NULL)
    {std::printf("Impossible d'ouvrir le fichier %s en ecriture \n",nomFichier);exit(1);}

  for(int i=0;i<nb;i++)
    {
      int temp;
      std::printf("Entier suivant : ");
      std::scanf("%d",&temp); //Ou utilisez la version robuste de saisie d'un int
      std::fwrite(&temp,sizeof(int),1,fd);
    }
  std::fclose(fd);
}

void litFichierEntiersBin(const char * nomFichier)
//preconditions : nomFichier chaine de caracteres designant le nom d'un fichier
//                contenant des entiers ecrits sous forme binaire
//postcondition : affichage du contenu du fichier sur la sortie standard
{
  FILE *fd=std::fopen(nomFichier,"r");
  if (fd==NULL)
    {std::printf("Impossible d'ouvrir le fichier %s en lecture \n",nomFichier);exit(1);}

  int temp;
  while(std::fread(&temp,sizeof(int),1,fd)==1)
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
  creeFichierEntiersBin(fich,n);
  std::printf("Affichage des entiers contenu dans le fichier\n");
  litFichierEntiersBin(fich);
  return 0;
}

#endif
