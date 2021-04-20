import smtplib
from email.header import Header
from email.mime.text import MIMEText
from email.utils import formataddr, parseaddr

def sendEmail(content, subject, sender, receivers, user, passwd, smtp_host, mail_format='plain'):
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
    

if __name__ == '__main__':
    # 第三方 SMTP 服务
    mail_host = "smtp.163.com"             # SMTP服务器
    mail_user = "yeting2938@163.com"     # 用户名
    mail_passwd = "yt151435411"            # 授权密码，非登录密码
    sender = 'yeting2938@163.com'
    receivers = ['yeting2938@qq.com', 'yeting2938@hust.edu.cn']  # 接收邮件，可设置为你的QQ邮箱或者其他邮箱
    content = '''
    <h4>单元格间距="0":</h4>
    <table border="1" cellspacing="0">
        <thead>
            <tr>
                <td> Titleeeeee </td>
                <td> Priceeeeee </td>
            </tr>
        </thead>
        
        <tbody>
            <tr>
                <td>AAA</td>
                <td>$53</td>
            </tr>
            <tr>
                <td>BBB</td>
                <td>$75</td>
            </tr>
        </tbody>
    </table>
    '''
    subject = '测试邮件'  # 邮件主题
    sendEmail(content, subject, sender, receivers, mail_user, mail_passwd, mail_host, 'html')
