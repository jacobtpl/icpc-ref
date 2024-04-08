alias c='g++ -Wall -Wconversion -Wfatal-errors -g -std=c++14 \
	-fsanitize=undefined,address'
xmodmap -e 'clear lock' -e 'keycode 66=less greater' #caps = <>
setxkbmap -option caps:escape # bind caps to esc
setxkbmap -option # remove options