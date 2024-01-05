import requests

def query_plain(text, url="http://bern2.korea.ac.kr/plain"):
    return requests.post(url, json={'text': text}).json()

if __name__ == '__main__':
    text = "a subset of clostridiaceae were depleted in stool corresponding with baseline adenomas and diabetes mellitus , while desulfovibrio was enriched both in stool and in mucosal biopsies ."
    print(query_plain(text))
    print(query_plain(text)['annotations'][0]['mention'])
    text.split(" ")