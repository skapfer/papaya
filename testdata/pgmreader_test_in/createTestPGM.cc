#include <fstream>
#include <cmath>

int main()
{
  std::ofstream outfile;
  outfile.open("test-16b-ascii.pgm");
  if (!outfile.good())
    return -1;
  outfile << "P2" << std::endl;
  outfile << "640 480" << std::endl;
  outfile << "65535" << std::endl;
  int mx=320; int my=240;
  for (int y=0; y<480; y++)
  for (int x=0; x<640; x++)
  {
    const double d=std::sqrt((x-mx)*(x-mx)+(y-my)*(y-my));
    outfile << (unsigned int)(65535*(std::max(1-d/mx,0.0))) << " ";
  }
  outfile.close();
  return 0;
}
