{-# LANGUAGE DoRec #-}
import Prelude hiding (catch)
import Types
import Eval
import Parse
import Control.Monad.Error
import Control.Exception (catch)
import Data.Map (toList)
import System.IO
import System.IO.Error hiding (catch)
import GHC.Conc (numCapabilities)

printWelcome :: IO ()
printWelcome = mapM_ putStrLn
  ["scheme522 interpreter version 0.707",
   "running on " ++ show numCapabilities ++ " cores",
   "type \":help\" for a list of toplevel commands",
   ""]

printHelp :: IO ()
printHelp = mapM_ putStrLn
  ["Available top-level commands:",
   "  :help       Shows this screen.",
   "  :bindings   Prints all the bindings in the current environment.",
   "              Bindings higher on the stack have higher priority.",
   "  :quit       Leaves the interpreter. You may also use Ctrl-C (Windows)",
   "              or Ctrl-D (*nix)\n"]

printPrefix :: IO ()
printPrefix = putStr "scheme522> "

printBindings :: Env -> IO ()
printBindings = mapM_ (\(a,b) -> putStrLn $ a ++ ": " ++ show b) . toList

printVal :: LispExpr -> IO ()
printVal expr = putStrLn ("=> " ++ show expr)

evalOne :: Env -> LispExpr -> LispMonad Env
evalOne env (List [Atom "define", Atom var, binding]) =
  do binding' <- eval binding env
     liftIO (printVal binding')
     return $ bind var binding' env

evalOne env (List [Atom "definerec", Atom var, body]) =
  do rec val <- eval body env'
         env' <- return $ bind var val env
     liftIO (printVal val)
     return env'

evalOne env (List [Atom a, List (Atom f : args), binding]) | a == "define"
                                                          || a == "definerec" =
  evalOne env desugared
  where desugared = List [Atom a,
                          Atom f,
                          List [Atom "lambda", List args, binding]]

evalOne env l@(List [Atom "load", expr]) =
  do path <- eval expr env
     case path of
       Primitive (String path) ->
         do exprs <- parseFile path
            evalMany env exprs
       _ -> throwError (TypeError l)

evalOne env expr =
  do val <- eval expr env
     liftIO (printVal val)
     return env

evalMany = foldM evalOne

runLine :: String -> Env -> LispMonad Env
runLine l env =
  do expr <- parseLine l
     evalMany env expr

main :: IO ()
main = do printWelcome
          loop emptyEnv
          where eofHandler e = if isEOFError e
                               then putStrLn "" >> return ":quit"
                               else ioError e
                loop env =
                  do printPrefix
                     hFlush stdout
                     line <- getLine `catch` eofHandler
                     case line of
                       ":help"     -> printHelp >> loop env
                       ":quit"     -> return ()
                       ":bindings" -> printBindings env >> loop env
                       _           ->
                         do result <- runLisp $ runLine line env
                            case result of
                              Left err -> print err >> loop env
                              Right v  -> loop v
