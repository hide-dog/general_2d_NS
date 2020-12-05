# general_2d_NS
コードの使い方  
1　general_2d_NSを好きなところに置く  
2　general_2d_NS.jl内のsrc_pathを変更する  
　　srcフォルダーの配下までを絶対パスで指定
3　pre.jlのセル数等を変更し，実行  
4　windowsのコマンドプロンプトで下記を実行しスレッド数を指定  
  ただし，格子数が少ないため，2か1でよい  
　　set JULIA_NUM_THREADS=2  
5　general_2d_NS.jlを実行  
6　post.jlを実行  
7　post_resultフォルダーに結果があるので，paraviewを使って可視化  
  
This software is released under the MIT License, see LICENSE.
