# See the status of your files (untracked, modified, staged)
git status
# Add all new/modified files (except ignored ones) to the staging area
git add .
# Commit the staged changes to your local repository with a descriptive message
git commit -m "Initial commit of ToCoSiAn project"
\end{tiny}
\item \textbf{Link to GitHub:}
\begin{enumerate}
\item Go to GitHub and create a new \textbf{private} repository (e.g., \texttt{ToCoSiAn-Project}).
\item GitHub will provide you with commands to link your local repository to this new remote one. They will look like this:
\end{enumerate}
\begin{minted}[
frame=lines,
framesep=2mm,
baselinestretch=1.2,
bgcolor=secondarycolor,
fontsize=\small,
breaklines=true
]{bash}
# Adds a remote named "origin" pointing to your GitHub URL
git remote add origin https://github.com/swarupdas2004/ToCoSiAn-Project.git
# Renames the default branch to "main" (modern convention)
git branch -M main
# Pushes your "main" branch to the remote "origin"
git push -u origin main