# Delta asc reader

JEOL RESONANCE社のNMR測定・解析ソフト「Delta」で作成できる「Generic ASCII」ファイル(*.asc + *.hdr)をnmrglueに読み込むためのライブラリです。

# 必要なライブラリと制限

* numpy
* scipy
* nmrglue>=0.8

Python 3.6以上で動作することを確認しています。現在、2Dまでのデータのみサポートしています。

# インストール方法
以下のようにすることでpipからインストール可能です;
```
$ pip install git+https://github.com/makki547/delta_asc_reader
```

## 使い方
ライブラリをインポートして`read`関数を使ってascファイルを読み込んでください。その際、ascファイルと同じディレクトリにhdrファイルが置かれている必要があります。
読み込まれたデータは`nmrglue.pipe.read`でNMRPipeファイルを読み込んだときと同様の形式、つまりヘッダ情報のディクショナリとnp.array型のデータを戻り値としています。

また、`get_array_values`関数を使うことで、アレイ測定した際のパラメータリストをリストとして取得できます。この際、SI接頭語(マイクロやミリ、キロなど)は削除されます。
例えば、アレイ測定でリスト`{0[ms], 1[ms], ... , 1000[ms]}`をy次元にセットした場合、`get_array_values('hoge.asc', 'y')`とすることで`[0, 0.001, 1.0]`というPythonリストを取得できます。

## 使用例

```py
import delta_asc_reader as delta
import nmrglue as ng

dic, data = delta.read('hoge.asc')
tau_intervals = delta.get_array_values('hoge.asc', 'y') #get y_list data

#Handle the data as same as dic, data from ng.pipe.read 
uc = ng.pipe.make_uc(dic, data, dim = 1)

```

## その他
このライブラリを使用して得られた結果、生じた損失等について、作者は一切の責任を負いません。自己責任で使用してください。
