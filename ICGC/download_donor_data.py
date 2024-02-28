import sys
import selenium
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.action_chains import ActionChains
import pandas as pd
import time
import os
import shutil

d = pd.read_csv(sys.argv[1], sep="\t")
fids = d["ICGC Donor"]
fids = list(fids)

base_url = "https://dcc.icgc.org/donors/"
path = "/home/user/Downloads"
o="/home/user/Downloads/icgc_download_data"

browser = webdriver.Firefox()

folder_names = []
total_files = len(fids)
missed = []

download_files = False
download_summary = True
sumd ={}

for x, i in enumerate(fids):
    print("On: ", x+1, "/", total_files, end="\r")
    url = base_url+i
    browser.get(url)
    wait = WebDriverWait(browser, 20)
    time.sleep(3)
    if download_files == True:
        # dl = browser.find_element(By.XPATH, '/html/body/div[3]/div[2]/span/article/div[2]/section/div[2]/div/button')
        dl = browser.find_element_by_css_selector(".t_button")
        wait.until(EC.visibility_of(dl))
        dl.click()
        time.sleep(1.5)
        #element = WebDriverWait(browser, 10).until_not(EC.presence_of_element_located((By.CSS_SELECTOR, ".modal-body > div:nth-child(2) > table:nth-child(2)")))
        #try:
        tab = browser.find_element_by_css_selector(".modal-body > div:nth-child(2) > table:nth-child(2)")
        wait.until(EC.visibility_of(tab))
        for row in tab.find_elements_by_xpath(".//tr"):
            #wait.until(EC.visibility_of(row))
            #coordinates = row.location_once_scrolled_into_view
            #browser.execute_script('window.scrollTo({}, {});'.format(coordinates['x'], coordinates['y']))
            #ActionChains(browser).move_to_element(row).perform()
            time.sleep(0.5)
            #row.click()
            browser.execute_script("arguments[0].click();", row)
        time.sleep(0.5)
        dlt = browser.find_element_by_xpath("/html/body/div[6]/div/div/div/div[3]/button[2]")
        dlt.click()
        time.sleep(5)
        files = (file for file in os.listdir(path) if os.path.isfile(os.path.join(path, file)))
        for file in files:
            os.rename(os.path.join(path, file), os.path.join(path, i+".tar"))
            time.sleep(0.1)
            shutil.move(os.path.join(path, i+".tar"), os.path.join(o, i+".tar"))
            time.sleep(0.1)
       # except:
       #     missed.append(i)
       #     print(missed)
       #     print(len(missed))
    if download_summary == True:
        sumd[i] = {}
        #summ = browser.find_element_by_xpath("/html/body/div[3]/div[2]/span/article/div[2]/section/div[1]/table")
        summ = browser.find_element_by_css_selector("div.half:nth-child(1) > table:nth-child(2)")
        wait.until(EC.visibility_of(summ))
        for row in summ.find_elements_by_xpath(".//tr"):
            sumd[i][row.find_element_by_xpath(".//th").text] = row.find_element_by_xpath(".//td").text
            #print(sumd)
        #print(sumd)

summary = pd.DataFrame(sumd)
summary.to_csv("summary.csv", index=False)
summary.to_csv("summary_no_quote.csv", index=False, quotechar=None)

