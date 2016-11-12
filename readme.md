#ICP
ICP
複数のレーザから得られた点群データ(Points Cloud)を繋ぐため  
点群間の最適な変異パラメータをICPで算出する  

ICPは以下を修正が改善する限り繰返す    
平行移動は各点群データの重心の差とする  
点群データを重心からの相対位置とする  
点群データ間の分散が最小になるSVD（特異値分解）で回転角度を推定  　 
修正が改善しない場合は終了する

### 点群データ　オリジナル◇点を回転・平行移動・ノイズ付与した＋点群にシフト　　
![image](https://cloud.githubusercontent.com/assets/20177544/20238332/e5fda37a-a92b-11e6-9c9c-a7aa8744213f.png)
### シフトした点群をICPによりオリジナル◇近傍に復元した点群□　かなり近接している  
![image](https://cloud.githubusercontent.com/assets/20177544/20238334/effdd91c-a92b-11e6-9878-59b128b96d27.png)
