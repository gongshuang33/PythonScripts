import requests 
import os
import random
import datetime
import json


BASEDIR = os.path.dirname(__file__)

# 日期
def getTimeStamp():
    return str(datetime.datetime.now().strftime("%Y%m%d.%H%M.%s"))

def getUrl(url):
    try:
        print(getTimeStamp())
        res = requests.get(url)
        templates_path = os.path.join(BASEDIR, "temps")
        if not os.path.exists:
            os.mkdir(templates_path)
        tempfile = os.path.join(templates_path, "respnse.html." + getTimeStamp())
        open(tempfile, 'w', encoding='utf-8').write(res.text)
        return {'status':0, "info":'', 'data':tempfile }
    except Exception as e:
        e = str(e) + str(e.__traceback__.tb_lineno)
        return {'status':-1, "info": str(e), 'data':'' }

if __name__ == "__main__":
    url = 'https://kg.qq.com/node/personal_v2?uid=639d95872c2d3e8b33&shareUid=&chain_share_id=nocy2BSG02WmKK1OVCGe7E7WSOpLJ-ZKGsqSVbGyoKU&pageId=homepage_guest'
    # print(json.dumps(getUrl(url)))
    print(getTimeStamp())