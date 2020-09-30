# create repository
git init
git remote add origin https://github.com/zhihuanli/Experimental-Method-in-Nuclear-Physics.git

#rutines for updates
git pull origin master  #将远程仓库同步到本地
git add [XXX]    #添加文件到暂存区    万能命令 git add --all
git commit -m "commit message"    #提交修改
git push -u origin master    #上传到网站 
