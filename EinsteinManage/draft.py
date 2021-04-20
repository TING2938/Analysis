#!/usr/bin/env python3

import os
import subprocess as sub
import numpy as np
import time
import smtplib
from email.header import Header
from email.mime.text import MIMEText

# for group
allGroup = {
    'fg': [
        'fg', 'yting', 'lzeng', 'zhaoc', 'sunchen', 'fjm', 'dqguo', 'myzhu', 'weiyan', 'wtao', 'wrx',
        'shuaiw', 'zzy', 'jkzhang', 'txy123', 'hljiang', 'bisheng', 'pfj', 'jxpeng', 'mtm', 'hhzhang', 'lz', 'zxwang'
    ],
    'ls': [
        'song', 'lw', 'hgq', 'ww', 'yxliu'
    ],
    'lh': [
        'luhao', 'shaosj', 'chenwz', 'lijinda', 'chenyang', 'zhangzhiyou'
    ],
    'pwu': [
        'pwu', 'suntao'
    ],
    'xbli': [
        'xbli35', 'ycong', 'xiongrui'
    ],
    'others': [
        'Guest', 'tes', 'dxh', 'guest', 'acm-user', 'lzh', 'xubo', 'chenlc', 'xky', 'yexiaoming', 'wukaidong', 'kyx', 'fhs', 'xxx', 'heshumo', 'jiangyankun', 'jack', 'icg_tarig', 'wangjun', 'tianjun', 'xkliu', 'dongxian', 'chukaixin'
    ]
}


mailList = {
    'fg': 'gfeng@hust.edu.cn',
    'ls': 'ls...'
}


# for mail setting
mail_host = 'smtp.163.com'
mail_user = 'hustitpeinstein@163.com'
mail_passwd = 'CHEIHNNEXCGYRORB'
sender = 'hustitpeinstein@163.com'
receivers = ["yeting2938@hust.edu.cn"]
# receivers = ["D201980386@hust.edu.cn", "mchen@hust.edu.cn", "yeting2938@hust.edu.cn"]

warningCount = {}
groupState = {}
groupStorage = {}
for key in allGroup.keys():
    warningCount[key] = 0
    groupState[key] = True
    groupStorage[key] = 0.0

def run(cmd: str, check: bool = False):
    h = sub.run(cmd, shell=True, stdout=sub.PIPE, stderr=sub.PIPE, check=check)
    return h.stdout.decode('utf-8')

def sendMail(content, subject, sender, receivers, user, passwd, smtp_host, mail_format='plain'):
    message = MIMEText(content, mail_format, 'utf-8')
    message['From'] = sender
    message['To'] = ",".join(receivers)
    message['Subject'] = subject
    try:
        smtpObj = smtplib.SMTP_SSL(smtp_host, 465)
        smtpObj.login(user, passwd)
        smtpObj.sendmail(sender, receivers, message.as_string())
    except smtplib.SMTPException as e:
        print(e)


while True:
    dfAvailable = float(run('df | egrep "/state/partition1$"').split()[4][0:-1])

    if dfAvailable > 85:
        dfcontent = ""
        for line in run("df -lh").split("\n")[1:]:
            dfcontent += "<tr>"
            for word in line.split():
                dfcontent += f'<td style="width: 15.5216%;">{word}</td>'
            dfcontent += "</tr>"
        content = open('EinsteinWarning.html.config', encoding='utf-8').read()
        content = content.replace("MAINTEXT_DF", dfcontent)
        content = content.replace("DataAndTime", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        sendMail(content, "Einsten Warning", sender, receivers, mail_user, mail_passwd, mail_host, "html")

        # os.system('cd /state/partition1/home && du -l --max-depth=1 > storage_for_all_user.dat')
        storage = np.loadtxt("/state/partition1/home/storage_for_all_user.dat", dtype=str)
        allMapStorage = {}
        for i in storage:
            allMapStorage[i[1][2:]] = float(i[0]) / 1e6  # GB

        for teacher, students in allGroup.items():
            groupStorage[teacher] = 0.0

            teacherTotalUsed = 0.0
            for student in students:
                teacherTotalUsed += allMapStorage[student]
                groupStorage[teacher] = teacherTotalUsed

            if teacherTotalUsed > 1.5 * 1000:  # 1.5TB
                warningCount[teacher] += 1
            else:
                warningCount[teacher] = 0

            if warningCount[teacher] >= 10:
                if groupState[teacher] == True:
                    os.system("echo disable permission")
                    os.system("echo send mail to teacher here")
                    groupState[teacher] = False
                else:
                    os.system("echo send mail to teacher here")

            if warningCount[teacher] == 0:
                if groupState[teacher] == False:
                    os.system("echo enable permission")
                    os.system("echo mail to teacher here")
                    groupState[teacher] = True

        totalMailContent = open('DetailedWarning.html.config', encoding='utf-8').read()
        content_df = ""
        content_group = ""
        content_student = ""

        for line in run("df -lh").split("\n")[1:]:
            content_df += "<tr>"
            for word in line.split():
                content_df += f'<td style="width: 15.5216%;">{word}</td>'
            content_df += "</tr>"

        for teacher, storage in groupStorage.items():
            content_group += f'''
            <tr>
                <td style="width: 15.5216%;">{teacher}</td>
                <td style="width: 15.5216%;">{storage:8.3f}</td>
            </tr>
            '''

        for teacher, students in allGroup.items():
            content_student += f'''
            <tr>
                <td style="width: 31.0432%; text-align: center;" colspan="2"><strong>(In group <span style="color: #e03e2d;">{teacher}</span>)</strong></td>
            </tr>
            '''
            for student in students:
                content_student += f'''
                <tr>
                    <td style="width: 15.5216%;">{student}</td>
                    <td style="width: 15.5216%;">{allMapStorage[student]:8.3f}</td>
                </tr>
                '''

        totalMailContent = totalMailContent.replace("MAINTEXT_DF", content_df)
        totalMailContent = totalMailContent.replace("MAINTEXT_Group", content_group)
        totalMailContent = totalMailContent.replace("MAINTEXT_Student", content_student)
        totalMailContent = totalMailContent.replace("DataAndTime", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

        sendMail(totalMailContent, "Einsten Detailed Warning", sender,
                 receivers, mail_user, mail_passwd, mail_host, "html")

    else:
        for teacher in allGroup.keys():
            warningCount[teacher] = 0
            if groupState[teacher] == False:
                os.system("echo enable permission here")
                os.system("echo mail to teacher here")
                groupState[teacher] = True

    time.sleep(2*3600)  # 2 hours
