# kosis_fetch.py - Downloads KOSIS tables through the KOSIS Open API, one JSON file per year. Usage: python3 kosis_fetch.py <orgId> <tblId> <first year> <last year> [OBJ_INDEX=item name ...]
#  Tables used: 350 DT_35006_N009 (long-term care certified persons by district, sex and grade), DT_35006_N036 (national certifications by sex and age),
#  DT_35006_N005 and DT_35006_N037 (applications), 117 DT_11761_N006/N009/N002/N007 (registered disability, screening comparison only).
#  The API key is read from the KOSIS_API_KEY environment variable. Output: ST_KOSIS_RAW/<tblId>/<year>.json
import os,json,urllib.request,urllib.parse,sys,time
KEY=os.environ['KOSIS_API_KEY']
org,tbl,y0,y1=sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4])
restr=dict(a.split('=',1) for a in sys.argv[5:])
OUT=os.path.join(os.environ['ST_KOSIS_RAW'],tbl); os.makedirs(OUT,exist_ok=True)
def get(url,q):
    for k in range(3):
        try: return json.loads(urllib.request.urlopen(url+urllib.parse.urlencode(q),timeout=300).read().decode('utf-8'))
        except Exception as e: err=e; time.sleep(3)
    return {'err':type(err).__name__}
m=get('https://kosis.kr/openapi/statisticsData.do?',dict(method='getMeta',type='ITM',apiKey=KEY,orgId=org,tblId=tbl,format='json',jsonVD='Y'))
objs=[];items=[]
for x in m:
    if x['OBJ_ID']=='ITEM': items.append(x['ITM_ID']); continue
    if x['OBJ_ID'] not in [o[0] for o in objs]: objs.append((x['OBJ_ID'],[]))
    dict(objs)[x['OBJ_ID']].append((x['ITM_ID'],x['ITM_NM']))
q0=dict(method='getList',apiKey=KEY,orgId=org,tblId=tbl,itmId='ALL',prdSe='Y',format='json',jsonVD='Y')
for i,(oid,lst) in enumerate(objs,1):
    if str(i) in restr:
        ids=[c for c,n in lst if n.strip() in restr[str(i)].split('|')]
        q0[f'objL{i}']='+'.join(ids)+'+' if len(ids)>1 else ids[0]
    else: q0[f'objL{i}']='ALL'
for y in range(y0,y1+1):
    f=f'{OUT}/{y}.json'
    if os.path.exists(f) and os.path.getsize(f)>100: continue
    q=dict(q0,startPrdDe=str(y),endPrdDe=str(y))
    d=get('https://kosis.kr/openapi/Param/statisticsParameterData.do?',q)
    if isinstance(d,dict): print(tbl,y,'ERR',str(d)[:120]); continue
    json.dump(d,open(f,'w'),ensure_ascii=False); print(tbl,y,len(d)); time.sleep(0.5)
