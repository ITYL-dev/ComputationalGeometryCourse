#ifdef COMPILE_CLASS_START

#include <iostream>
#include <vector>
#include <map>

class LaClasse
{
 public : // fonctions membres publiques
  //Construction, conversion, affectation et destruction
  LaClasse() : l(0)
  {std::cout << "LaClasse::LaClasse()\n";}
  LaClasse(const LaClasse & lc) : l(lc.l)
  {std::cout << "LaClasse::LaClasse(const LaClasse & )\n";}
  LaClasse(int i) : l(i) 
  {std::cout << "LaClasse::LaClasse(int)\n";}
  operator bool() const
  {std::cout << "LaClasse::operator bool() const\n"; return true;}
  ~LaClasse()
  {std::cout << "~LaClasse::LaClasse()\n";}
  const LaClasse & operator=(const LaClasse & c)
  {l=c.l; std::cout << "LaClasse::operator=(const LaClasse &)\n"; return *this;}
  //Autre fonction membre
  LaClasse F(LaClasse);
  // Declaration fonction exterieure amie
  friend std::ostream & operator << (std::ostream & os, const LaClasse & lc);
private : // donnees membres privees
  int l;
};

LaClasse F(LaClasse vv)
{
  std::cout << " in F \n";
  return 8;
}

LaClasse LaClasse::F(LaClasse v)
{
  std::cout << " in LaClasse::F \n";
  return ::F(v);
}

std::ostream & operator << (std::ostream & os, const LaClasse & lc)
{
  os << " in ostream << LaClasse "<< lc.l << std::endl;
  return os;
}

// Testez et analysez la séquence d'appels aux fonctions membres 
// de LaClasse dans le programme suivant :


int main()
{
  std::cout << "Les classes en C++ (constructeurs, destructeurs, affectations) \n";
  LaClasse c1;
    LaClasse cc1(c1); // meme chose que LaClasse cc1=c1;
    LaClasse cv1(5);  // meme chose que LaClasse cv1=5;
  std::cout << std::endl;
  c1=cv1;
  std::cout << std::endl;
    
  std::cout << "Les conversions, les fonctions, les fonctions membres et les valeurs de retour \n";
  cc1=8;
  std::cout << std::endl;
  if(cv1)
    {
      cv1=F(10);
      cv1=F(c1);
      cv1=c1.F(cc1);
    }
  std::cout << std::endl;
    
  std::cout << "Tableaux \n";
  LaClasse tab[3];
  tab[1]=cv1;
  LaClasse *pc=new LaClasse(tab[0]);
    *pc=c1; //Meme chose que pc[0]=c1;
  delete pc;
  std::cout << std::endl;
    
  std::cout << "Les vectors (tableaux dynamiques) \n";
  std::vector<int> vv;
  int i;
  std::cout << "Entrez des entiers (0 pour finir):\n";
  do {
    std::cin >> i;
    vv.push_back(i);
    } while (i);
  std::cout << "vv contient " << vv.size() << " entiers.\n";
  std::cout << "Dans la case " << vv.size()/2 <<" il y a " << vv[vv.size()/2] << "\n";
  std::cout << std::endl;
    
  std::cout << "Les maps (tables associatives ou dictionnaires) \n";
  std::map<char,int> mam;
  mam.insert( std::pair<char,int>('c',5) );
  mam.insert( std::pair<char,int>('e',10) );
    std::cout << "En fait insert retourne un iterateur et un boolen pour indiquer ou l'element a ete insere (booleen a true) ou bien ou il etait deja (booleen a false) \n";
  std::pair<std::map<char,int>::iterator,bool> ret;
  ret = mam.insert( std::pair<char,int>('n',20) );
  if (ret.second==false) {
    std::cout << "L'element y etait deja, ";
    std::cout << "avec la valeur " << ret.first->second << '\n';
  }
    
  std::map<char,int>::iterator it;
  std::cout << "mam contient:\n";
  for (it=mam.begin(); it!=mam.end(); ++it)
    std::cout << it->first << " => " << it->second << '\n';
  return 0; 
}

#endif