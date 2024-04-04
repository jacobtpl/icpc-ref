set ts=2 sw=2 ai cin nu rnu udf udir=~/.vim/udir 
sy on
set cul ru nowrap wmnu sc is bs=indent,eol,start cino=g0
" Select region and then type :Hash to hash your selection.
" Useful for verifying that there aren't mistypes.
ca Hash w !cpp -dD -P -fpreprocessed \| tr -d '[:space:]' \
 \| md5sum \| cut -c-6

nmap <F8> :w <Bar> !g++ -std=c++20 -DLOCAL %<CR>
nmap <F9> :w <Bar> !g++ -std=c++20 -DLOCAL % && ./a.out<CR>

autocmd FileType python set sw=4 ts=4 sts=4 et nocin si
autocmd FileType python nmap <F8> :w <Bar> !python3 "%"<CR>
autocmd FileType python nmap <F9> :w <Bar> !python3 -i "%"<CR>

" To map caps lock -> escape --
" setxkbmap -option caps:escape
