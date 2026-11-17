#include "../utilities/template.h"

#include "../../content/number-theory/ModInverse.h"

int main() {
	initInv();
	rep(i,1,LIM) assert(inv[i] * i % mod == 1);
	cout<<"Tests pass!"<<endl;
}
