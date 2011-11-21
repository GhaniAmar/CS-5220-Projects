module Types where

import Prelude hiding (catch)
import qualified Data.Map as M
import Data.List (intercalate)
import Control.Monad.Error
import Control.Exception
import Text.ParserCombinators.Parsec.Error (ParseError)

-- Type declarations -----------------------------------------------------------

type Symbol = String

data Primitive = Integer Integer
               | Double Double
               | String String
               | Bool Bool
                 deriving (Eq, Ord)

data LispExpr = Atom Symbol
              | List [LispExpr]
              | Dotted [LispExpr] LispExpr
              | Primitive Primitive
              | Closure LispExpr Symbol Env
                deriving (Eq, Ord)

type Env = M.Map Symbol LispExpr

emptyEnv :: Env
emptyEnv = M.empty

data LispErr = TypeError LispExpr
             | IOError IOError
             | SyntaxError ParseError
             | NumericError LispExpr
             | Malformed LispExpr
             | UnboundValue LispExpr
             | NonFunctional LispExpr [LispExpr]
             | EmptyList LispExpr
             | OtherError String

type LispMonad = ErrorT LispErr IO

-- Type instances --------------------------------------------------------------

instance Show Primitive where
    show (Integer n) = show n
    show (Double n)  = show n
    show (String s)  = "\"" ++ s ++ "\""
    show (Bool b)    = if b then "#t" else "#f"

instance Show LispExpr where
    show (Atom s) = s
    show (List l) = "(" ++ intercalate " "  (map show l) ++ ")"
    show (Dotted l t) = "(" ++ intercalate " " (map show l) ++ " . " ++ show t ++ ")"
    show (Primitive p) = show p
    show (Closure _ _ _) = "<fun>"

instance Show LispErr where
  show (TypeError e)       = "Type error in expression " ++ show e
  show (IOError e)         = "IO error: " ++ show e
  show (SyntaxError e)     = "Syntax error: " ++ show e
  show (NumericError e)    = "Numeric error encountered in expression " ++ show e
  show (Malformed e)       = "Malformed expression: " ++ show e
  show (UnboundValue e)    = "Unbound value: " ++ show e
  show (NonFunctional e _) = "Non functional value applied: " ++ show e
  show (EmptyList e)       = "Illegal list operation on empty list: " ++ show e
  show (OtherError s)      = "Error: " ++ s

instance Error LispErr where
  strMsg = OtherError

-- Functions on these types ----------------------------------------------------

lookupEnv :: Symbol -> Env -> LispMonad LispExpr
lookupEnv "#t" _ = return . Primitive . Bool $ True
lookupEnv "#f" _ = return . Primitive . Bool $ False
lookupEnv s env = case M.lookup s env of
  Nothing -> throwError . UnboundValue $ Atom s
  Just v  -> return v

bind :: Symbol -> LispExpr -> Env -> Env
bind = M.insert

runLisp :: LispMonad a -> IO (Either LispErr a)
runLisp m = catch (runErrorT m) (return . Left . IOError)

malformed :: LispExpr -> LispMonad LispExpr
malformed = throwError . Malformed

toBool :: LispExpr -> LispMonad Bool
toBool (Primitive (Bool b)) = return b
toBool e = throwError (TypeError e)