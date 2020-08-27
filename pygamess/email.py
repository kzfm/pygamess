from email.mime.text import MIMEText
import smtplib, os, ssl
#import traceback
context = ssl.create_default_context()

def smtplib_email(email_body, receivers, subject, smtpconfig, debug_level=0):
    #If exception is None, this function is deemed to be reporting a succesful run, else it is deemed to
    #be reporting an error
    use_ssl = smtpconfig.get("ssl")
    if use_ssl:
        smtp = smtplib.SMTP_SSL(smtpconfig['server'], smtpconfig['port'], context=context)
    else:
        smtp = smtplib.SMTP(smtpconfig['server'], smtpconfig['port'])
    with smtp:
       #smtp.set_debuglevel(1)
        smtp.set_debuglevel(debug_level)
        if use_ssl:
            smtp.login(smtpconfig["username"], smtpconfig["password"])
            #smtp.starttls(context=context)
        if isinstance(receivers, str):
            receivers = (receivers, )#Create a singleton tuple, tuple because tuples are more efficient than lists
        msg = MIMEText(email_body)
        msg['From'] = smtpconfig['username']
        msg['To'] = ','.join(receivers)
        msg['Subject'] = subject
        smtp.send_message(msg)
